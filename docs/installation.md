## COWBAT Installation

### Dependencies

* Linux system
* [Conda](https://conda.io/docs/user-guide/install/linux.html)
* [Docker](https://www.docker.com/)

#### Conda method

The way I install conda:

```
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
bash miniconda.sh -b -p $HOME/miniconda
export PATH="$HOME/miniconda/bin:$PATH"
conda config --set always_yes yes
conda update -q conda
```

The easiest way to install COWBAT is to download the source code [GitHub Link](https://github.com/OLC-Bioinformatics/COWBAT.git)

```
git clone https://github.com/OLC-Bioinformatics/COWBAT.git
cd COWBAT
export PATH="/path/to/repository/COWBAT:$PATH"
conda env create -f environment.yml
source activate cowbat
```

#### Docker method

Docker must already be installed

The docker image relies on conda to install all the dependencies, so the cowbat environment must be sourced within 
the container prior to launch. The supplied command below launches container, and immediately sources the environment, and runs the 
pipeline, but it is also possible to run those commands separately from within the container. For additional details on the run
command, please see [the tutorial](tutorial.md).

```
git clone https://github.com/OLC-Bioinformatics/COWBAT.git
cd COWBAT
docker build -t cowbat:latest .
docker run -it --name cowbat --rm cowbat:latest /bin/bash -c "source activate cowbat && assembly_pipeline.py -s /path/to/sequences -r /path/to/database"
```

### Databases

Use the database setup script included in OLCTools. This will download and set up all required databases.

NOTE: If you want rMLST databases, you must contact Keith Jolley (keith.jolley@zoo.ox.ac.uk)
for an account, and for the necessary keys.

```
python -m databasesetup.database_setup -d /PATH/TO/DESIRED/LOCATION -c /PATH/TO/RMLST/CREDENTIALS
```

### Testing

[Unit tests](tests.md)

