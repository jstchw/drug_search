const BASE_URL = 'https://api.fda.gov/drug/event.json'


const generatePath = (searchTerm, searchOptions, countType) => {
    const fromDate = `20040101`
    const toDate = formatDate(new Date())
    const whatToCount = {
        reaction: 'patient.reaction.reactionmeddrapt.exact',
        receiveDate: 'receivedate',
        drug: 'patient.drug.activesubstance.activesubstancename.exact'
    }

    let searchParts = [`(receivedate:[${fromDate}+TO+${toDate}])`]

    if (searchTerm) {
        // Join search terms with '+' and add to searchParts
        searchParts.push(`${searchOptions.searchBy}:${searchTerm.join('+')}`);
    }

    if (searchOptions.sex) {
        searchParts.push(searchOptions.sex)
    }

    if (searchOptions.age && searchOptions.age[0] && searchOptions.age[1]) {
        searchParts.push(`patient.patientonsetage:[${searchOptions.age[0]}+TO+${searchOptions.age[1]}]`)
    }

    if (searchOptions.country) {
        searchParts.push(`occurcountry:"${searchOptions.country}"`)
    }

    return `${BASE_URL}?search=${searchParts.join('+AND+')}&count=${whatToCount[countType]}`
}



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

export const getEventsFromDrugs = async (searchTerm, searchOptions) => {
    const url = generatePath(searchTerm, searchOptions, 'reaction')
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

export const getEventsOverTime = async (searchTerm, searchOptions) => {
    const url = generatePath(searchTerm, searchOptions, 'receiveDate')
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

export const getDrugsFromEvents = async (searchTerm, searchOptions) => {
    const url = generatePath(searchTerm, searchOptions, 'drug')
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