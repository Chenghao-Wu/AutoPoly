# CLI

The `autopoly` command is installed with the package (`autopoly=AutoPoly.cli:main`). All output is JSON to stdout, errors to stderr; exit code 0 on success, 1 on failure.

## Commands

```bash
autopoly info                     # all options, limits, and example configs
autopoly validate <config.json>   # validate a config file
autopoly validate --stdin         # validate from stdin
autopoly generate <config.json>   # generate a polymer system
autopoly generate --stdin         # generate from stdin
autopoly describe <smiles>        # describe a SMILES string
```

## Examples

```bash
# Discover what AutoPoly can do
autopoly info

# Validate, then generate
autopoly validate pe.json && autopoly generate pe.json

# Pipe a config directly
echo '{"type":"atomistic","name":"pe","polymers":[{"chain_num":2,"sequence":["CC[*]","[*]CC[*]","[*]CC"]}]}' \
  | autopoly generate --stdin
```

See the [Agent API & CLI guide](../guides/agent.md) for the full config schema.

## Module

::: AutoPoly.cli
