# Automated superfamily renaming
## What and why?
This repo exists to automate the renaming of current superfamilies in CATH according to [InterPro guidelines](./InterPro_guidelines.md).
***

### Automation

Most of the job is done by the custom python module [CATH-parser](CATH_parser.py).

#### [Renamed superfamilies](./results/renamed_superfamilies.tsv)
Has the easiest cases, which refer to the first part of the guidelines. Fixes errors such as trailing stops, semicolons, lowercase start (except mRNA and such) etc. Whenever a fix happens a one-letter comment is added.

##### Comment codes

* **R** - removed excessive capital letters
* **S** - replace semicolon
* **L** - lowercase start
* **T** - trailing stop
* **C** - other stop(except between digits)
* **I** - ii to II conversion, for example photosystem II
* **O** - capitalising first letter in various bacteria *(only coccus so far)*
* **N** - capitalising single letters, mostly n/c-terminal and cytochrome C

![replace](./plots/replacement.png)
 ***
#### [Flagged superamilies]('./results/flagged,tsv')

Superfamilies that have something wrong with the name but are too complicated to fix automatically.

* **P** - punctuation (currently underscore, plus and colon)
* **F** - forbidden words
* **S** - bad start
* **W** - unprreferred words

![flag](./plots/flags.png)
***
#### [Duplicated superfamilies](./results/duplicates.tsv)
Has a list of superfamilies with identical names. Mostly it is pairs however, there are a lot of superfamilies *named helix hairpin bin* and *single helix bin*.

### Manual curation

Even though automation deals with the vast majority of errors, there are still a number of issues that need to be dealt with manually. This is why the [manual curation list](./manual_curation_flags.tsv) exists. It contains suggestions for renaming that are easier done manually rather than via automation.
