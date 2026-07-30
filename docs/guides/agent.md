# Agent API & CLI

AutoPoly exposes a **config-driven JSON interface** designed for scripts, pipelines, and AI agents: plain dicts in, JSON-serializable results out, no exceptions raised. The same interface backs the `autopoly` command-line tool and optional LangChain tool wrappers.

## The five functions

| Function | Purpose |
|---|---|
| `info()` | Discover all options: force fields, topologies, bead-spring choices, limits, and complete example configs. **Call this first.** |
| `validate(config)` | Check a config without writing anything. Returns errors, warnings, and fix suggestions. |
| `generate(config)` | Validate, then build the system. Returns output paths and metadata. |
| `describe_smiles(smiles)` | Atom counts, elements, molecular weight, and wildcard info for a SMILES string. |
| `suggest_force_field(smiles_list)` | Suggest compatible force fields from the elements present. |

```python
from AutoPoly import agent

# 1. Discover what's available
options = agent.info()

# 2. Validate a config (never writes to disk)
result = agent.validate(config)
if not result.success:
    print(result.errors, result.suggestions)

# 3. Generate
result = agent.generate(config)
if result.success:
    print(result.data_file, result.files_created)
```

## Config schema

Every config has a `type`: `"atomistic"` or `"bead_spring"`.

### Atomistic

```json
{
  "type": "atomistic",
  "name": "pe_system",
  "output_dir": "./output",
  "force_field": "oplsaa",
  "placement_method": "mc_random",
  "use_mc_chain_growth": true,
  "polymers": [
    {
      "chain_num": 10,
      "sequence": ["CC[*]", "[*]CC[*]", "[*]CC"],
      "topology": "linear",
      "tacticity": "atactic"
    }
  ],
  "molecules": [
    {"count": 100, "smiles": "O", "name": "water"}
  ]
}
```

Polymer sequences use [complement SMILES](complement-smiles.md); molecule SMILES must **not** contain wildcards. `force_field` accepts the six values in the [Force Fields guide](force-fields.md); `placement_method` is `"grid"` or `"mc_random"` (see [MC Placement](mc-placement.md)).

### Bead-spring

```json
{
  "type": "bead_spring",
  "name": "melt",
  "output_dir": "./output",
  "n_chains": 100,
  "bead_types": [{"name": "A"}, {"name": "B", "epsilon": 1.2}],
  "sequence": [["A", 50], ["B", 50]],
  "topology": "linear",
  "bond_style": "harmonic",
  "pair_style": "lj",
  "generation_method": "saw"
}
```

Fields map one-to-one onto `BeadSpringPolymer` — see the [Bead-Spring guide](bead-spring.md).

### Results

`generate()` returns a `GenerationResult` with `success`, `output_dir`, `data_file`, `files_created`, and `metadata` (force field, chain/bead counts, box size). On failure it carries `errors` and `warnings` — it never raises. `validate()` returns a `ValidationResult` with `success`, `errors`, `warnings`, and `suggestions` (including "did you mean" hints for misspelled options).

## Command line

The `autopoly` command ships with the package. All output is JSON on stdout, errors on stderr; exit code 0 on success, 1 on failure.

```bash
autopoly info                        # all options, limits, example configs
autopoly validate config.json        # check a config file
autopoly validate --stdin            # ... or pipe it
autopoly generate config.json        # build the system
autopoly describe "CC[*]"            # inspect a SMILES string
```

```bash
echo '{"type":"bead_spring","name":"t","n_chains":5,"bead_types":[{"name":"A"}],"sequence":[["A",20]]}' \
  | autopoly validate --stdin
```

## LangChain tools

With the optional extra installed (`pip install -e ".[agent]"`), AutoPoly provides ready-made LangChain tools with Pydantic schemas:

```python
from AutoPoly.tools import get_autopoly_tools

tools = get_autopoly_tools()
# autopoly_info, autopoly_generate_atomistic,
# autopoly_generate_bead_spring, autopoly_describe_smiles
```

These plug into any LangChain/DeepAgents agent loop — the intended pattern is: the agent calls `autopoly_info` to learn the config format, drafts a config, validates it, then generates.

## Resource limits

`agent.info()` reports the hard limits (also enforced by `validate()`): max DOP 10 000, max sequence length 10 000, max 100 unique monomers per sequence.

## See also

- [agent API reference](../reference/agent.md) — function signatures from source
- [tools API reference](../reference/tools.md) — LangChain wrappers
- [CLI reference](../reference/cli.md) — the `autopoly` command
