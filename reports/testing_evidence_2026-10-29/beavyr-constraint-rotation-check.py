import json
import math

def sub(a, b): return [x - y for x, y in zip(a, b)]
def dot(a, b): return sum(x * y for x, y in zip(a, b))
def cross(a, b): return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]
def norm(v): return math.sqrt(dot(v, v))
def distance(a, b): return norm(sub(a, b))
def dihedral(p):
    b0, b1, b2 = sub(p[0], p[1]), sub(p[2], p[1]), sub(p[3], p[2])
    axis = [x / norm(b1) for x in b1]
    v = [x - dot(b0, axis)*y for x, y in zip(b0, axis)]
    w = [x - dot(b2, axis)*y for x, y in zip(b2, axis)]
    return math.degrees(math.atan2(dot(cross(axis, v), w), dot(v, w)))

# The same four carbon coordinates used by the failing production test.
before = [[-1.92, 0.26, 0.0], [-0.64, -0.57, 0.0], [0.64, 0.26, 0.0], [1.92, -0.57, 0.0]]
axis = sub(before[2], before[1])
axis = [x/norm(axis) for x in axis]
v = sub(before[3], before[2])
phi = math.pi / 3
perpendicular = cross(axis, v)
rotated = [v[i]*math.cos(phi) + perpendicular[i]*math.sin(phi) + axis[i]*dot(axis, v)*(1-math.cos(phi)) for i in range(3)]
after = before[:3] + [[before[2][i]+rotated[i] for i in range(3)]]
changes = [abs(distance(before[a], before[b]) - distance(after[a], after[b])) for a, b in [(0, 1), (1, 2)]]
angle_change = abs((dihedral(after) - dihedral(before) + 180) % 360 - 180)
assert max(changes) < 1e-12
assert abs(angle_change - 60) < 1e-10
print(json.dumps({'fixed_bond_length_changes_angstrom': changes, 'dihedral_change_degrees': angle_change,
                  'conclusion': 'These bond-length constraints allow a 60-degree torsion change; they do not make butane rigid.'}, indent=2))
