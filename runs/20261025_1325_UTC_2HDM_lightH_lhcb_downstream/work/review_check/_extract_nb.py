import json
p = r"E:/AAAwork/LLP-AGI/references/tools/FORESEE/Models/2HDM_H/2HDM_H.ipynb"
nb = json.load(open(p, encoding="utf-8"))
for i, c in enumerate(nb["cells"]):
    src = "".join(c.get("source", []))
    if "9.68544e-13" in src:
        print("=== cell", i, "===")
        print(src)
