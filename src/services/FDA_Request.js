const BASE_URL = 'https://api.fda.gov/drug/event.json'

const formatDate = (date) => {
    const year = date.getFullYear()
    const month = date.getMonth() + 1
    const day = date.getDate()
    return `${year}${month < 10 ? `0${month}` : month}${day < 10 ? `0${day}` : day}`
}

export const getDrugEventsSearch = async (searchTerm, searchType) => {
    const fromDate = `20040101`
    const toDate = formatDate(new Date())
    const count = 'patient.reaction.reactionmeddrapt.exact'

    const url = `${BASE_URL}?search=(receivedate:[${fromDate}+TO+${toDate}])+AND+${searchType}:"${encodeURIComponent(searchTerm)}"&count=${count}`
    const response = await fetch(url)
    return await response.json()
}

export const getDrugEventsCount = async (searchType, searchTerm, limit) => {
    const url = `${BASE_URL}?count=${searchType}:${searchTerm}&limit=${limit}`
    const response = await fetch(url)
    return await response.json()
}