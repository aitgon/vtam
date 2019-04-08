import tempfile



tempdir = tempfile.mkdtemp()

# Old tax_assign parameters
order = [100.0, 97.0, 95.0, 90.0]
# order = [100.0, 97.0, 95.0, 90.0, 85.0, 80.0]
taxonomic_levels = {"family": 5, "order": 4, "genus": 3, "species": 2, "subspecies": 1}

# New tax_assign parameters
rank_hierarchy =['no rank', 'phylum', 'superclass', 'class', 'subclass', 'infraclass', 'superorder', 'order',
                 'suborder', 'infraorder', 'family', 'subfamily', 'genus', 'subgenus', 'species', 'subspecies']


