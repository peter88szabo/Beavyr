use super::*;
use std::time::{SystemTime, UNIX_EPOCH};

fn directory(label: &str) -> PathBuf {
    let path = std::env::temp_dir().join(format!(
        "beavyr python {label} {} {}",
        std::process::id(),
        SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    std::fs::create_dir_all(&path).unwrap();
    path
}

#[test]
fn an_explicit_interpreter_preserves_paths_and_arguments_without_a_shell() {
    let settings = Settings {
        interpreter: "/some other disk/my venv/bin/python".into(),
        ..Default::default()
    };
    let command = settings.command(["a script.py", "$(never execute this)"]);
    assert_eq!(command.get_program(), "/some other disk/my venv/bin/python");
    assert_eq!(
        command.get_args().collect::<Vec<_>>(),
        ["a script.py", "$(never execute this)"]
    );
    assert_eq!(command.get_envs().count(), 0);
}

#[test]
fn backend_settings_are_saved_separately_and_round_trip() {
    let root = directory("settings");
    let pyscf = Settings {
        activation: "source '/external disk/pyscf/bin/activate'".into(),
        ..Default::default()
    };
    let psi4 = Settings {
        activation: "micromamba activate psi4env".into(),
        ..Default::default()
    };
    save_to(&root.join("pyscf.ron"), &pyscf).unwrap();
    save_to(&root.join("psi4.ron"), &psi4).unwrap();
    assert_eq!(load_from(&root.join("pyscf.ron")).unwrap(), pyscf);
    assert_eq!(load_from(&root.join("psi4.ron")).unwrap(), psi4);
    std::fs::write(root.join("broken.ron"), "invalid settings").unwrap();
    assert!(load_from(&root.join("broken.ron")).is_err());
    std::fs::remove_dir_all(root).unwrap();
}

#[test]
fn powershell_commands_quote_interpreter_and_arguments_as_literals() {
    let settings = Settings {
        interpreter: "C:\\User's env\\python.exe".into(),
        shell: "powershell.exe".into(),
        activation: "micromamba activate psi4env".into(),
    };
    let command = settings.command(["a script.py", "$unexpanded"]);
    let args: Vec<_> = command
        .get_args()
        .map(|a| a.to_string_lossy().into_owned())
        .collect();
    assert_eq!(&args[..3], ["-NoLogo", "-NonInteractive", "-Command"]);
    assert!(args[3].contains("micromamba activate psi4env\n"));
    assert!(args[3].contains("& 'C:\\User''s env\\python.exe' 'a script.py' '$unexpanded'"));
}

#[cfg(unix)]
#[test]
fn activation_failure_stops_before_python_and_is_not_reported_ready() {
    let settings = Settings {
        shell: "/bin/bash".into(),
        activation: "false".into(),
        ..Default::default()
    };
    let root = directory("failure");
    let output = settings
        .command(["-c", "print('CALCULATION_STARTED')"])
        .env("HOME", &root)
        .output()
        .unwrap();
    assert!(!output.status.success());
    assert!(!String::from_utf8_lossy(&output.stdout).contains("CALCULATION_STARTED"));
    assert!(String::from_utf8_lossy(&output.stderr).contains("activation failed"));
    assert!(!super::super::probe_with_settings("json", &settings).is_ready());
    std::fs::remove_dir_all(root).unwrap();
}

#[cfg(unix)]
#[test]
fn activation_loads_shell_hooks_and_changes_only_the_child_environment() {
    let root = directory("shell hooks");
    std::fs::write(
        root.join(".bashrc"),
        "micromamba() { [ \"$1\" = activate ] && export BEAVYR_TEST_ACTIVE=\"$2\"; }\n",
    )
    .unwrap();
    let settings = Settings {
        shell: "/bin/bash".into(),
        activation: "micromamba activate 'an env anywhere'".into(),
        ..Default::default()
    };
    let parent_value = std::env::var_os("BEAVYR_TEST_ACTIVE");
    let output = settings
        .command(["-c", "import os; print(os.environ['BEAVYR_TEST_ACTIVE'])"])
        .env("HOME", &root)
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        String::from_utf8_lossy(&output.stdout).trim(),
        "an env anywhere"
    );
    assert_eq!(std::env::var_os("BEAVYR_TEST_ACTIVE"), parent_value);
    std::fs::remove_dir_all(root).unwrap();
}

#[cfg(unix)]
#[test]
fn a_venv_in_an_arbitrary_path_is_used_for_both_probe_and_execution() {
    let root = directory("venv ' with spaces");
    let prefix = root.join("environment");
    let status = Command::new(super::super::PYTHON)
        .args(["-m", "venv", "--without-pip"])
        .arg(&prefix)
        .status()
        .unwrap();
    assert!(status.success());
    let activation_path = prefix
        .join("bin/activate")
        .to_string_lossy()
        .replace('\'', "'\\''");
    let settings = Settings {
        activation: format!("source '{activation_path}'"),
        shell: "/bin/bash".into(),
        ..Default::default()
    };
    let result = super::super::probe_with_settings("json", &settings);
    let super::super::PythonEnv::Ready {
        executable,
        prefix: found,
        ..
    } = result
    else {
        panic!("{result:?}")
    };
    assert_eq!(Path::new(&found), prefix);
    let output = settings
        .command(["-c", "import sys; print(sys.executable)"])
        .output()
        .unwrap();
    assert!(output.status.success());
    assert_eq!(String::from_utf8_lossy(&output.stdout).trim(), executable);
    let direct = Settings {
        interpreter: prefix.join("bin/python").display().to_string(),
        ..Default::default()
    };
    assert!(super::super::probe_with_settings("json", &direct).is_ready());
    std::fs::remove_dir_all(root).unwrap();
}

#[cfg(unix)]
#[test]
#[ignore = "Uses the explicitly provided BEAVYR_TEST_PYTHON_ACTIVATION command"]
fn real_psi4_environment_can_be_activated_and_run() {
    let settings = Settings {
        activation: std::env::var("BEAVYR_TEST_PYTHON_ACTIVATION").unwrap(),
        ..Default::default()
    };
    let environment = super::super::probe_with_settings("psi4", &settings);
    assert!(environment.is_ready(), "{environment:?}");
    let scratch = directory("psi4 calculation");
    let output = settings.command(["-c", "import psi4; psi4.set_num_threads(1); psi4.core.be_quiet(); mol = psi4.geometry('0 1\\nH 0 0 0\\nH 0 0 0.74'); print('BEAVYR_ENERGY', psi4.energy('hf/sto-3g', molecule=mol))"]).current_dir(&scratch).output().unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(String::from_utf8_lossy(&output.stdout).contains("BEAVYR_ENERGY -1."));
    std::fs::remove_dir_all(scratch).unwrap();
}
