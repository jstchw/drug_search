const BASE_URL = 'https://api.fda.gov/drug/event.json'

const formatDate = (date) => {
    const year = date.getFullYear()
    const month = date.getMonth() + 1
    const day = date.getDate()
    return `${year}${month < 10 ? `0${month}` : month}${day < 10 ? `0${day}` : day}`
}

const processDrugEvents = (data) => {
    const termCountDict = {}

    const totalCount = data.results.reduce((total, item) => {
        termCountDict[item.term] = item.count
        return total + item.count
    }, 0)
    return { totalCount, termCountDict }
}

export const getDrugEventsSearch = async (searchTerm, searchType) => {
    const fromDate = `20040101`
    const toDate = formatDate(new Date())
    const count = 'patient.reaction.reactionmeddrapt.exact'

    const url = `${BASE_URL}?search=(receivedate:[${fromDate}+TO+${toDate}])+AND+${searchType}:"${encodeURIComponent(searchTerm)}"&count=${count}`
    console.log(url)
    try {
        const response = await fetch(url)
        if (!response.ok) {
            throw new Error(`Error: ${response.status}: ${response.statusText}`)
        }
        const data = await response.json()
        return {result: processDrugEvents(data)}
    } catch (error) {
        return {error: error.message}
    }
}

export const getEventsOverTime = async (searchTerm, searchType) => {
    const fromDate = `20040101`
    const toDate = formatDate(new Date())
    const count = 'receivedate'
    const url = `${BASE_URL}?search=(receivedate:[${fromDate}+TO+${toDate}])+AND+${searchType}:"${encodeURIComponent(searchTerm)}"&count=${count}`
    try {
        const response = await fetch(url)
        if (!response.ok) {
            throw new Error(`Error: ${response.status}: ${response.statusText}`)
        }
        const data = await response.json()
        return {result: data}
    } catch (error) {
        return {error: error.message}
    }
}