from pathlib import Path
DATA={'expected.json': '{\n  "objective": "(x-1)^2+(y-2)^2",\n  "equality": "x+y=3",\n  "initial_point": [\n    0.5,\n    2.5\n  ],\n  "bounds": [\n    -5,\n    5\n  ],\n  "optimal_x": [\n    1,\n    2\n  ],\n  "optimal_objective": 0,\n  "solution_tolerance": 0.001,\n  "constraint_tolerance": 1e-05\n}', 'problem.txt': '0.5 2.5 -5 5\n'}
root=Path(__file__).resolve().parents[1]/"examples"
root.mkdir(exist_ok=True)
for name,text in DATA.items(): (root/name).write_text(text,encoding="utf-8")
print("Generated analytical problems in",root)
