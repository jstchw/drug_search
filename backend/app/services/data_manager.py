import json


class DataManager:
    _instance = None

    @classmethod
    def get_instance(cls):
        if cls._instance is None:
            cls._instance = cls()
        return cls._instance

    def __init__(self):
        if DataManager._instance is not None:
            raise Exception("This class is a singleton!")


        # Data initialization for active substance suggestion
        with open('data/suggestions/substances.json', 'r') as file:
            self.substances = json.load(file)

        # Data initialization for brand name suggestion
        with open('data/suggestions/products.json', 'r') as file:
            self.products = json.load(file)

        with open('data/suggestions/side_effects.json', 'r') as file:
            self.side_effects = json.load(file)

        # Data initialization for DRUGBANK JSON
        with open('data/drugbank/drugs.json', 'r') as file:
            self.drug_json = json.load(file)

        with open('data/pubmed_article_data/pubmed_data.json', 'r') as file:
            self.pubmed_data = json.load(file)
