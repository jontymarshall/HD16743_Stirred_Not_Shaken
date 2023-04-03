This readme file contains instructions that will be useful in the github repository, or can be read as stand a
lone -- instructions on how to run the entire analysis, either through Makefile, runall.py or runall.R

```
.
├── analysis            <- Directory to work on your analysis
├── data
│   ├── external	<- Data from third party sources.
│   ├── intermediate	<- Intermediate data that has been transformed.
│   ├── processed	<- The final, canonical data sets for modeling.
│   └── raw		<- The original, immutable data dump.
├── docs		<- Space for documentation
├── .env		<- Store your secrets and config variables in a special file
├── env			<- Will contain the Python executable files and installed libraries for your virtualenv environment
├── .gitignore		<- Avoids uploading data, credentials, outputs, system files etc
├── LICENCE
├── models		<- Trained and serialized models, model predictions, or model summaries
├── paper		<- Manuscript or other dissemination files
│   ├── biblio.bib
│   ├── definitions.tex
│   ├── figures		<- Generated graphics and figures to be used in reporting
│   ├── paper_template.org
│   └── section_or_appendix
│       └── description
│           ├── figures
│           └── table.tex
├── README.md		<- The top-level README for developers using this project.
├── references		<- Data dictionaries, manuals, etc.
├── requirements.txt	<- Install the environment dependencies with: `pip install -r requirements.txt`
├── runall.py		<- execute the pipeline in python
└── tox.ini		<- Automate testing, cf. https://tox.readthedocs.io/en/latest/.
```

