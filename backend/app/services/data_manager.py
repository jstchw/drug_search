import json
from app.utils import load_data


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

        # Data initialization for SUBSTANCES AND PRODUCTS
        self.substances = load_data('data/suggestions/suggestion_db.csv', ',')
        self.products = load_data('data/suggestions/products_fda.csv', ',')

        # Data initialization for DRUGBANK JSON
        with open('data/json/drugs_cleaned_pretty_without_nulls.json', 'r') as file:
            self.drug_json = json.load(file)
