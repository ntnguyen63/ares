# ares
Polishing pipeline using racon, medaka, and pilon. 

## Install
```
git clone https://github.com/Jearce/ares.git
cd ares
mamba env create --file ares-env.yaml 
chmod +x ares.py
conda activate ares-env
```
Lastly, add ares directory to your $PATH to make ares.py globally available. For example in your .bashrc file:
```
export PATH="$PATH:/path/to/ares"
```

## Help
```
ares.py --help
```

```
Usage: ares.py [OPTIONS] DRAFTS NANOPORE R1 R2 BUSCO_LINEAGE T

Arguments:
  DRAFTS         [required]
  NANOPORE       [required]
  R1             [required]
  R2             [required]
  BUSCO_LINEAGE  [required]
  T              [required]

Options:
  --install-completion  Install completion for the current shell.
  --show-completion     Show completion for the current shell, to copy it or
                        customize the installation.

  --help                Show this message and exit.
```
