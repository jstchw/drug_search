age_groups = {
    'infant': {
        'min': 0,
        'max': 2,
        'fda_key': 2,
    },
    'child': {
        'min': 2,
        'max': 12,
        'fda_key': 3,
    },
    'adolescent': {
        'min': 12,
        'max': 18,
        'fda_key': 4,
    },
    'adult': {
        'min': 18,
        'max': 65,
        'fda_key': 5,
    },
    'elderly': {
        'min': 65,
        'max': 200,
        'fda_key': 6,
    }
}

sex_groups = {
    'unknown': {
        'fda_key': 0,
    },
    'male': {
        'fda_key': 1,
    },
    'female': {
        'fda_key': 2,
    }
}

count_properties_from_search_type = {
    'generic_name': [
        'effect',
        'treatment_disorder',
    ],
    'brand_name': [
        'effect',
        'treatment_disorder',
    ],
    'side_effect': 'generic_name',
}

search_fields: dict = {
    'generic_name': ['drug'],
    'brand_name': ['drug'],
    'side_effect': [
        'effect',
        'treatment_disorder',
    ]
}