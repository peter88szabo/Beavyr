//! Streaming GIF89a with per-frame median-cut palettes and bounded memory.
use std::io::{self, Write};

pub struct Gif<W: Write> {
    writer: W,
    width: u16,
    height: u16,
    fps: u16,
    frames: u64,
}
impl<W: Write> Gif<W> {
    pub fn new(mut writer: W, width: u16, height: u16, fps: u16, repeat: bool) -> io::Result<Self> {
        if width == 0 || height == 0 || fps == 0 || fps > 100 {
            return Err(io::Error::other("Invalid GIF dimensions or frame rate"));
        }
        writer.write_all(b"GIF89a")?;
        writer.write_all(&width.to_le_bytes())?;
        writer.write_all(&height.to_le_bytes())?;
        writer.write_all(&[0x70, 0, 0])?;
        if repeat {
            writer.write_all(b"\x21\xff\x0bNETSCAPE2.0\x03\x01\x00\x00\x00")?;
        }
        Ok(Self {
            writer,
            width,
            height,
            fps,
            frames: 0,
        })
    }
    pub fn frame(&mut self, rgba: &[u8]) -> io::Result<()> {
        if rgba.len() != self.width as usize * self.height as usize * 4 {
            return Err(io::Error::other("Wrong GIF frame size"));
        }
        let (palette, pixels) = quantize(rgba);
        self.writer.write_all(&[0x21, 0xf9, 4, 4])?;
        let delay = (((self.frames + 1) * 100 / self.fps as u64)
            - (self.frames * 100 / self.fps as u64)) as u16;
        self.frames += 1;
        self.writer.write_all(&delay.to_le_bytes())?;
        self.writer.write_all(&[0, 0])?;
        self.writer.write_all(&[0x2c, 0, 0, 0, 0])?;
        self.writer.write_all(&self.width.to_le_bytes())?;
        self.writer.write_all(&self.height.to_le_bytes())?;
        self.writer.write_all(&[0x87])?;
        self.writer.write_all(&palette)?;
        self.writer.write_all(&[8])?;
        let compressed = weezl::encode::Encoder::new(weezl::BitOrder::Lsb, 8)
            .encode(&pixels)
            .map_err(|e| io::Error::other(format!("GIF compression: {e:?}")))?;
        for block in compressed.chunks(255) {
            self.writer.write_all(&[block.len() as u8])?;
            self.writer.write_all(block)?;
        }
        self.writer.write_all(&[0])
    }
    pub fn finish(mut self) -> io::Result<()> {
        self.writer.write_all(&[0x3b])?;
        self.writer.flush()
    }
}
#[derive(Clone)]
struct Color {
    rgb: [u8; 3],
    count: u64,
    key: usize,
}
fn key(rgb: &[u8]) -> usize {
    ((rgb[0] as usize >> 3) << 10) | ((rgb[1] as usize >> 3) << 5) | (rgb[2] as usize >> 3)
}
fn range(colors: &[Color]) -> (u8, usize) {
    (0..3)
        .map(|axis| {
            (
                colors.iter().map(|c| c.rgb[axis]).max().unwrap_or(0)
                    - colors.iter().map(|c| c.rgb[axis]).min().unwrap_or(0),
                axis,
            )
        })
        .max()
        .unwrap()
}
fn quantize(rgba: &[u8]) -> ([u8; 768], Vec<u8>) {
    let mut bins = vec![[0u64; 4]; 32768];
    for p in rgba.chunks_exact(4) {
        let b = &mut bins[key(p)];
        for c in 0..3 {
            b[c] += p[c] as u64;
        }
        b[3] += 1;
    }
    let colors: Vec<_> = bins
        .iter()
        .enumerate()
        .filter(|(_, b)| b[3] > 0)
        .map(|(key, b)| Color {
            rgb: [
                (b[0] / b[3]) as u8,
                (b[1] / b[3]) as u8,
                (b[2] / b[3]) as u8,
            ],
            count: b[3],
            key,
        })
        .collect();
    let mut boxes = vec![colors];
    while boxes.len() < 256 {
        let Some(index) = (0..boxes.len())
            .filter(|&i| boxes[i].len() > 1)
            .max_by_key(|&i| range(&boxes[i]).0)
        else {
            break;
        };
        let mut colors = boxes.swap_remove(index);
        let axis = range(&colors).1;
        colors.sort_by_key(|c| c.rgb[axis]);
        let total: u64 = colors.iter().map(|c| c.count).sum();
        let mut sum = 0;
        let mut cut = 1;
        for (i, c) in colors.iter().enumerate().take(colors.len() - 1) {
            sum += c.count;
            cut = i + 1;
            if sum >= total / 2 {
                break;
            }
        }
        let second = colors.split_off(cut);
        boxes.push(colors);
        boxes.push(second);
    }
    let mut palette = [0u8; 768];
    let mut lookup = vec![0u8; 32768];
    for (i, colors) in boxes.iter().enumerate() {
        let total = colors.iter().map(|c| c.count).sum::<u64>().max(1);
        for axis in 0..3 {
            palette[3 * i + axis] = (colors
                .iter()
                .map(|c| c.rgb[axis] as u64 * c.count)
                .sum::<u64>()
                / total) as u8;
        }
        for c in colors {
            lookup[c.key] = i as u8;
        }
    }
    (
        palette,
        rgba.chunks_exact(4).map(|p| lookup[key(p)]).collect(),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn fractional_centisecond_frame_rate_preserves_total_duration() {
        let mut bytes = vec![];
        let mut gif = Gif::new(&mut bytes, 1, 1, 24, false).unwrap();
        for _ in 0..24 {
            gif.frame(&[0, 0, 0, 255]).unwrap();
        }
        gif.finish().unwrap();
        let delays: Vec<_> = bytes
            .windows(8)
            .filter(|w| w[..4] == [0x21, 0xf9, 4, 4])
            .map(|w| u16::from_le_bytes([w[4], w[5]]))
            .collect();
        assert_eq!(delays.len(), 24);
        assert_eq!(delays.iter().sum::<u16>(), 100);
    }
    #[test]
    fn gif_stream_has_header_frames_repeat_delay_and_trailer() {
        let mut bytes = vec![];
        let mut gif = Gif::new(&mut bytes, 2, 1, 20, true).unwrap();
        gif.frame(&[209, 114, 4, 255, 235, 161, 0, 255]).unwrap();
        gif.frame(&[0, 0, 0, 255, 255, 255, 255, 255]).unwrap();
        gif.finish().unwrap();
        assert_eq!(&bytes[..6], b"GIF89a");
        assert_eq!(bytes.last(), Some(&0x3b));
        assert!(bytes.windows(11).any(|w| w == b"NETSCAPE2.0"));
        assert_eq!(
            bytes
                .windows(4)
                .filter(|w| *w == [0x21, 0xf9, 4, 4])
                .count(),
            2
        );
    }
    #[test]
    fn palette_preserves_distinct_element_colors_and_lzw_round_trips() {
        let rgba = [209, 114, 4, 255, 235, 161, 0, 255];
        let (palette, pixels) = quantize(&rgba);
        for (i, p) in rgba.chunks_exact(4).enumerate() {
            assert_eq!(
                &palette[pixels[i] as usize * 3..pixels[i] as usize * 3 + 3],
                &p[..3]
            );
        }
        let encoded = weezl::encode::Encoder::new(weezl::BitOrder::Lsb, 8)
            .encode(&pixels)
            .unwrap();
        assert_eq!(
            weezl::decode::Decoder::new(weezl::BitOrder::Lsb, 8)
                .decode(&encoded)
                .unwrap(),
            pixels
        );
    }
}
