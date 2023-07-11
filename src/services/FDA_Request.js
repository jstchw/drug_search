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
        searchTerm.map(term => searchParts.push(`${searchOptions.searchBy}:"${encodeURIComponent(term)}"`))
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
        const data = await response.json()
        return {result: processDrugEvents(data)}
    } catch (error) {
        throw new Error(error.message)
    }
}

export const getEventsOverTime = async (searchTerm, searchOptions) => {
    const url = generatePath(searchTerm, searchOptions, 'receiveDate')
    try {
        const response = await fetch(url)
        const data = await response.json()
        return {result: data}
    } catch (error) {
        throw new Error(error.message)
    }
}

export const getDrugsFromEvents = async (searchTerm, searchOptions) => {
    const url = generatePath(searchTerm, searchOptions, 'drug')
    try {
        const response = await fetch(url)
        const data = await response.json()
        return {result: processDrugEvents(data)}
    } catch (error) {
        throw new Error(error.message)
    }
}

const demographicDataCache = {}

export const getSideEffectsForDemographics = async (drug) => {
    const demographicGroups = [
        { name: 'Young Boys', sex: '1', age: [0, 18], def: 'Males aged 0-18 years'},
        { name: 'Young Girls', sex: '2', age: [0, 18], def: 'Females aged 0-18 years'},
        { name: 'Men', sex: '1', age: [19, 59], def: 'Males aged 19-59 years'},
        { name: 'Women', sex: '2', age: [19, 59], def: 'Females aged 19-59 years'},
        { name: 'Elderly Women', sex: '2', age: [60, 120], def: 'Females aged 60-120 years'},
        { name: 'Elderly Men', sex: '1', age: [60, 120], def: 'Males aged 60-120 years'}
    ];

    // If the data for this drug is already in the cache, return it
    if (demographicDataCache[drug]) {
        return demographicDataCache[drug];
    }

    const results = {};

    for (const group of demographicGroups) {
        const searchOptions = {
            searchBy: 'patient.drug.medicinalproduct',
            sex: `patient.patientsex:${group.sex}`,
            age: group.age,
        };

        const url = generatePath([drug], searchOptions, 'reaction')
        const response = await fetch(url)
        const data = await response.json()
        console.log('works')
        results[group.name] = {
            ...processDrugEvents(data),
            def: group.def,
            age: group.age
        }
    }

    // Save the results to the cache
    demographicDataCache[drug] = results

    return results;
};
