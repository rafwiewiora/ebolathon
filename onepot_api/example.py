"""Minimal onepot direct-API examples. Run: python example.py

Reads the key from ONEPOT_API_KEY. All three calls were verified live 2026-07-18.
"""
import os
from onepot import Client

client = Client(api_key=os.environ["ONEPOT_API_KEY"])

# 1) Sample a property-filtered, diverse slice of makeable CORE space (Scenario B).
space = client.sample_space(
    count=5,
    strategy="diverse",
    properties={"molecular_weight": {"min": 300, "max": 400}, "qed": {"min": 0.5}},
    include_properties=True,
)
print("seed:", space["seed"])
for m in space["molecules"]:
    print(f"  ${m['price_usd']:>3} {m['smiles']}")

# 2) Analog search with retrosynthetic decomposition.
resp = client.search(
    smiles_list=["c1ccc(NC(=O)c2ccccc2)cc1"],
    max_results=10,
    decompose=True,
)
print("credits_used:", resp["credits_used"])
q = resp["queries"][0]
for h in q["results"][:5]:
    print(f"  sim={h['similarity']:.3f} ${h['price_usd']} rxn={h.get('reaction_class', '-')} {h['smiles']}")

# 3) Synthon-window search: reuse a reaction_class from a prior decompose response,
#    tighten one position (exploit) and open another (explore).
#    filters = [{"reaction_class": "rxn_a46303e0", "bb_index": 0, "min_similarity": 0.8},
#               {"reaction_class": "rxn_a46303e0", "bb_index": 1, "max_similarity": 0.6}]
#    client.search(smiles_list=[anchor], bb_filters=filters, decompose=True)
