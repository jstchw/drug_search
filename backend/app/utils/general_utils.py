from constants import age_groups, sex_groups

def generate_param_list(params: dict) -> list:
    """
    Generate a list of parameters for an FDA search based on the given parameters.

    Args:
        params (dict): A dictionary containing search parameters such as 'search_type'

    Returns:
        list: A list of parameters for an FDA search.
    """
    param_list = []
    group_type = params.get('group_type')
    terms = params.get('terms')
    country = params.get('country')
    search_mode = params.get('search_mode')

    if group_type == 'sex':
        sex_param = params.get('sex')
        for _, age_range in age_groups.items():
            param = {
                'terms': terms,
                'sex': sex_param,
                'age': {
                    'min': age_range['min'],
                    'max': age_range['max']
                },
                'country': country,
                'search_mode': search_mode
            }
            param_list.append(param)
    elif group_type == 'age':
        age_param = params.get('age')
        for sex_group, _ in sex_groups.items():
            if sex_group != 'unknown':
                param = {
                    'terms': terms,
                    'sex': sex_group,
                    'age': age_param,
                    'country': country,
                    'search_mode': search_mode
                }
                param_list.append(param)

    return param_list