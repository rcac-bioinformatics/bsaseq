# Bioconda Recipe for bsaseq

## Submission instructions

1. Fork bioconda-recipes: https://github.com/bioconda/bioconda-recipes
2. Create branch: `git checkout -b bsaseq`
3. Copy this recipe to `recipes/bsaseq/meta.yaml`
4. Update SHA256 hash from PyPI:
   ```bash
   # Get SHA256 from PyPI
   curl -sL https://pypi.io/packages/source/b/bsaseq/bsaseq-1.0.0.tar.gz | sha256sum
   ```
5. Submit pull request

## Testing locally

```bash
# Install bioconda-utils
conda install -c bioconda bioconda-utils

# Lint recipe
bioconda-utils lint recipes/bsaseq

# Build recipe
bioconda-utils build recipes/bsaseq

# Test installation
conda create -n test-bsaseq -c ./output bsaseq
conda activate test-bsaseq
bsaseq --help
```

## After merge

Package will be available via:
```bash
conda install -c bioconda bsaseq
```

## Version updates

When releasing a new version:
1. Update `version` in meta.yaml
2. Update SHA256 hash
3. Reset build number to 0
4. Submit PR to bioconda-recipes

## Dependencies

The recipe uses `noarch: python` since bsaseq is pure Python.
The main binary dependency (htslib via cyvcf2) is handled by the
cyvcf2 conda package.
