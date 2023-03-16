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

/*
export const getDrugEventsSearch = async (searchType, searchTerm, limit) => {
    const url = `${BASE_URL}?search=${searchType}:${searchTerm}&limit=${limit}`
    const response = await fetch(url)
    return await response.json()
}

This is my function in FDA_Request.js. How do I structure it properly to retrieve data?

https://api.fda.gov/drug/event.json?search=(receivedate:[20040101+TO+20230314])+AND+patient.drug.activesubstance.activesubstancename:"modafinil"&limit=1

This is the API link. patient.drug.activesubstance.activesubstancename is a searchType, modafinil is a searchTerm and "1" is a limit
*/