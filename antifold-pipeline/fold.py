import json, os, subprocess, sys

F = os.path.expanduser("~/.claude/skills/nebius-scientific-ai/forge.py")

def a3m(path):
    d = json.load(open(path))
    return d["result"]["alignments"]["pdb70_220313"]["a3m"]["alignment"]

VH_MSA   = a3m("results/msa-vh/result.json")
VL_MSA   = a3m("results/msa-vl/result.json")
HER2_MSA = a3m("results/msa-her2/result.json")

VL   = "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK"
HER2 = open("her2.txt").read().strip()

def chains(vh):
    def poly(i, s, m):
        return {"id": i, "molecule_type": "protein", "sequence": s,
                "msa": {"msa_search": {"a3m": {"alignment": m, "format": "a3m"}}}}
    return [poly("H", vh, VH_MSA), poly("L", VL, VL_MSA), poly("T", HER2, HER2_MSA)]

designs = []
for block in open("designs/designs.fasta").read().split(">")[1:]:
    lines = block.strip().split("\n")
    designs.append((lines[0], "".join(lines[1:])))

os.makedirs("folds", exist_ok=True)
for name, vh in designs:
    args = {"polymers": chains(vh), "recycling_steps": 3, "diffusion_samples": 1}
    p = f"folds/{name}.json"
    json.dump(args, open(p, "w"))
    print("submitting", name)
    subprocess.run(["python3", F, "run", "boltz2", "--args-file", p,
                    "--label", f"fold-{name}", "--no-wait"])