const BASE_URL = 'https://api.fda.gov/drug/event.json'

export const getDrugEventsSearch = async (searchType, searchTerm, limit) => {
    const url = `${BASE_URL}?search=${searchType}:${searchTerm}&limit=${limit}`
    const response = await fetch(url)
    return await response.json()
}

export const getDrugEventsCount = async (searchType, searchTerm, limit) => {
    const url = `${BASE_URL}?count=${searchType}:${searchTerm}&limit=${limit}`
    const response = await fetch(url)
    return await response.json()
}