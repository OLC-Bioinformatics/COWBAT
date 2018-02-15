## COWBAT Installation

### Dependencies

* Linux system
* [Conda](https://conda.io/docs/user-guide/install/linux.html)

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

### Databases

Use the database_setup.py script included in the repository. This will download and set up all required databases.

NOTE: If you want rMLST databases, you must contact Keith Jolley (keith.jolley@zoo.ox.ac.uk)
for an account, and for the necessary keys.

```
python database_setup.py -d /PATH/TO/DESIRED/LOCATION 
```

### Testing

[Unit tests](tests.md)

