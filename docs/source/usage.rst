# Usage

## Installation

To install NSForest with pip: 

```shell
pip install nsforest
```

To install NSForest via Github: 

```shell
git clone https://github.com/JCVenterInstitute/NSForest.git
cd NSForest
conda env create -f nsforest.yml
conda activate nsforest
```

## Running

Running NSForest

```shell
python3 nsforest -a ${input_folder}/arguments_${prefix}.csv 
```

Parallelizing NSForest

```shell
python3 nsforest -a ${input_folder}/arguments_${prefix}.csv -c "${cluster}"
```

Creating recipes
----------------

To retrieve a list of random ingredients,
you can use the ``lumache.get_random_ingredients()`` function:

.. autofunction:: lumache.get_random_ingredients

The ``kind`` parameter should be either ``"meat"``, ``"fish"``,
or ``"veggies"``. Otherwise, :py:func:`lumache.get_random_ingredients`
will raise an exception.

.. autoexception:: lumache.InvalidKindError

For example:

>>> import lumache
>>> lumache.get_random_ingredients()
['shells', 'gorgonzola', 'parsley']

