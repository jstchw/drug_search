from typing import Optional
from datetime import date
from constants import fda_event_base_url, fda_search_sex, fda_search_types, fda_count_types, age_groups, sex_groups
import asyncio
from aiohttp import ClientSession, ClientConnectorError, ClientOSError, ServerDisconnectedError
from utils.general_utils import generate_param_list

async def fetch(session, url):
    try:
        async with session.get(url, timeout=10) as response:
            response.raise_for_status()
            return await response.json()
    except (ClientConnectorError, ClientOSError, ServerDisconnectedError):
        # Handle connection errors
        pass
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

def generate_fda_url(params: dict, entry_count: bool = False, limit = 50) -> str:
    """
    Generate an FDA URL based on the given parameters.

    Args:
        params (dict): A dictionary containing search parameters such as 'search_type
        entry_count (bool): If True, returns a URL to get the count of entries

    Returns:
        str: The generated FDA URL.
    """
    from_date = "20040101"
    current_date = date.today()
    to_date = date(current_date.year - 1, 12, 31).strftime("%Y%m%d")


    default_min_age, default_max_age = 0, 120
    min_age, max_age = default_min_age, default_max_age

    age_params = params.get('age')
    if age_params is not None:
        min_age = age_params.get('min')
        max_age = age_params.get('max')
    
        if min_age is None:
            min_age = default_min_age
        if max_age is None:
            max_age = default_max_age
    
    age_range = f"[{min_age}+TO+{max_age}]"

    search_parts = []
    sex_param = params.get('sex')
    if sex_param is not None:
        sex = fda_search_sex.get(sex_param)
        if sex:
            search_parts.append(sex)

    terms_params = params.get('terms')
    if terms_params is not None:
        search_parts.extend(
            f"({fda_search_types.get(term_obj['type'], '')}:{'+'.join(term_obj['term'].split())})"
            for term_obj in terms_params
            if term_obj['type'] in fda_search_types
        )
        what_to_count = fda_count_types.get(terms_params[0]['type'], '')

    country_param = params.get('country')
    if country_param is not None:
        search_parts.append(f"occurcountry:\"{country_param}\"")

    search_parts.append(f"patient.patientonsetage:{age_range}" if age_params is not None else '')
    search_parts.insert(0, f"(receivedate:[{from_date}+TO+{to_date}])")

    search_query = '+AND+'.join(part for part in search_parts if part)

    if entry_count:
        return f"{fda_event_base_url}?search={search_query}&limit=1"
    else:
        return f"{fda_event_base_url}?search={search_query}&count={what_to_count}&limit={limit}"
    

async def search_fda(params: dict, limit = 50) -> Optional[list]:
    """
    Search the FDA database for drug data based on the given parameters.

    Args:
        params (dict): A dictionary containing search parameters such as 'search_type'

    Returns:
        Optional[list]: A list of matched drug data, or None if an error occurred.
    """
    if params.get('view', '') == 'simple':
        result_url = generate_fda_url(params, limit=limit)

        async with ClientSession() as session:
            result_obj = await fetch(session, result_url)

            if result_obj is not None:
                # Return both results and total_count together
                return get_simple_breakdown_fda(result_obj['results'])
            
    elif params.get('view', '') == 'advanced':
        param_list = generate_param_list(params)
        url_param_tuples = [(generate_fda_url(param, limit=limit), param) for param in param_list]

        async with ClientSession() as session:
            result_objs = await asyncio.gather(
                *[fetch(session, url) for url, _ in url_param_tuples]
            )

            results_data = []
            for result_obj, (_, param) in zip(result_objs, url_param_tuples):
                if result_obj:
                    results_data.append({
                        'results': result_obj['results'],
                        'param': param
                    })

        return get_advanced_breakdown_fda(results_data, params)

    return None


def get_simple_breakdown_fda(data: list) -> dict:
    """
    Get a simple breakdown of the given FDA data.

    Args:
        data (list): A list of FDA data.
        total_count (int): The total count of matched entries.

    Returns:
        dict: A simple breakdown of the given FDA data.
    """
    
    categories = [entry['term'].capitalize() for entry in data]
    series_data = [entry['count'] for entry in data]

    series = [{
        'data': series_data,
        'name': 'Count'
    }]

    return {
        'categories': categories,
        'series': series,
        'total_count': sum(series_data),
    }


def get_advanced_breakdown_fda(data: list, params: dict) -> dict:
    """
    Get an advanced breakdown of the given FDA data.

    Args:
        data (list): A list of FDA data.
        params (dict): A dictionary of parameters used to generate the data.

    Returns:
        dict: An advanced breakdown of the given FDA data.
    """
    categories = []
    series_data = {}
    all_terms = set()

    if params.get('group_type', '') == 'sex':
        # Create categories based on age groups
        categories = [f"{age_group['min']} - {age_group['max']}" for age_group in age_groups.values()]
    elif params.get('group_type', '') == 'age':
        # Create categories based on sex groups (excluding 'unknown')
        categories = [sex_group for sex_group in sex_groups.keys() if sex_group != 'unknown']

    # Initialize series data with an empty dictionary for each category
    for category in categories:
        series_data[category] = {}

    # Iterate over each data entry
    for entry in data:
        category = None
        if params.get('group_type', '') == 'sex':
            # Get the age range from the param of the entry
            age_min = entry['param']['age']['min']
            age_max = entry['param']['age']['max']
            category = next((f"{age_group['min']} - {age_group['max']}" for age_group in age_groups.values()
                             if age_min == age_group['min'] and age_max == age_group['max']), None)
        elif params.get('group_type', '') == 'age':
            # Get the sex from the param of the entry
            sex = entry['param']['sex']
            category = sex if sex in categories else None

        if category is not None:
            for result in entry['results']:
                term = result['term']
                all_terms.add(term)
                series_data[category][term] = series_data[category].get(term, 0) + result['count']

    # Create the final series data
    final_series_data = []
    for term in all_terms:
        term_data = {
            'name': term.capitalize(),
            'data': [series_data[category].get(term, 0) for category in categories]
        }
        final_series_data.append(term_data)

    return {
        'categories': [category.capitalize() for category in categories],
        'series': final_series_data,
        'total_count': sum(sum(data) for data in zip(*[series['data'] for series in final_series_data])),
    }