age_groups = {
    'infant': {
        'min': 0,
        'max': 2,
        'fda_key': 2,
    },
    'child': {
        'min': 3,
        'max': 11,
        'fda_key': 3,
    },
    'adolescent': {
        'min': 12,
        'max': 17,
        'fda_key': 4,
    },
    'adult': {
        'min': 18,
        'max': 64,
        'fda_key': 5,
    },
    'elderly': {
        'min': 65,
        'max': 120,
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
        'effect'
    ],
    'brand_name': [
        'effect'
    ],
    'side_effect': 'generic_name',
}

search_fields: dict = {
    'generic_name': ['drug'],
    'brand_name': ['drug'],
    'side_effect': [
        'effect',
    ]
}

fda_event_base_url = 'https://api.fda.gov/drug/event.json'

fda_search_types = {
    'generic_name': 'patient.drug.openfda.generic_name',
    'brand_name': 'patient.drug.openfda.brand_name',
    'side_effect': 'patient.reaction.reactionmeddrapt',
}

fda_count_types = {
    'generic_name': 'patient.reaction.reactionmeddrapt.exact',
    'brand_name': 'patient.reaction.reactionmeddrapt.exact',
    'side_effect': 'patient.drug.activesubstance.activesubstancename.exact',
}

fda_search_sex = {
    'unknown': 'patient.patientsex:0',
    'male': 'patient.patientsex:1',
    'female': 'patient.patientsex:2',
}

#https://api.fda.gov/drug/event.json?
#search=(receivedate:[20040101+TO+20231231])
#+AND+(patient.reaction.reactionmeddrapt:fatigue)
#+AND+patient.patientonsetage:[0+TO+120]&
#count=patient.drug.activesubstance.activesubstancename.exact&limit=50