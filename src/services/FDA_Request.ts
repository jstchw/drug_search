import {DemographicGroups, FDARawData, SearchTypesMap, SearchOptions} from "../types";

const BASE_URL = 'https://api.fda.gov/drug/event.json'



const generatePath = (searchTerms: string[], searchOptions: SearchOptions, countType: keyof SearchTypesMap) => {
    const fromDate = `20040101`
    const toDate = formatDate(new Date())
    const whatToCount: SearchTypesMap = {
        reaction: 'patient.reaction.reactionmeddrapt.exact',
        receiveDate: 'receivedate',
        drug: 'patient.drug.activesubstance.activesubstancename.exact'
    }

    const searchParts = [`(receivedate:[${fromDate}+TO+${toDate}])`]

    // After the Search page has been translated to TS, work on the searchBy in searchOptions
    // It needs to be passed as a whole object, not just the value and then the value needs to be extracted
    if (searchTerms) {
        const encodedTerms = searchTerms.flat().join('+AND+')
        searchParts.push(`(${searchOptions.searchBy.value}:${encodedTerms})`)
    }

    if (searchOptions.sex) {
        searchParts.push(searchOptions.sex.value)
    }

    if (searchOptions.age && searchOptions.age.min.value && searchOptions.age.max.value) {
        searchParts.push(`patient.patientonsetage:[${searchOptions.age.min.value}+TO+${searchOptions.age.max.value}]`)
    }

    if (searchOptions.country) {
        searchParts.push(`occurcountry:"${searchOptions.country.value}"`)
    }

    return `${BASE_URL}?search=${searchParts.join('+AND+')}&count=${whatToCount[countType]}`
}



const formatDate = (date: Date) => {
    const year = date.getFullYear() - 1
    const month = date.getMonth() + 1
    const day = date.getDate()
    return `${year}${month < 10 ? `0${month}` : month}${day < 10 ? `0${day}` : day}`
}

const processDrugEvents = (data: FDARawData) => {
    const termCountDict: {[key: string]: number} = {}

    const totalCount = data.results.reduce((total, item) => {
        termCountDict[item.term] = item.count
        return total + item.count
    }, 0)
    return { totalCount, termCountDict }
}

export const getEventsFromDrugs = async (searchTerm: string, searchOptions: SearchOptions) => {
    const url = generatePath([searchTerm], searchOptions, 'reaction')
    try {
        const response = await fetch(url)
        const data = await response.json()
        return {result: processDrugEvents(data)}
    } catch (e) {
        if (e instanceof Error) {
            return {error: e.message}
        } else if (typeof e === 'string') {
            return {error: e}
        } else {
            return {error: 'Unknown error when searching by active substance or brand name'}
        }
    }
}

export const getEventsOverTime = async (searchTerm: string, searchOptions: SearchOptions) => {
    const url = generatePath([searchTerm], searchOptions, 'receiveDate')
    try {
        const response = await fetch(url)
        const data = await response.json()
        return {result: data}
    } catch (e) {
        if (e instanceof Error) {
            return {error: e.message}
        } else if (typeof e === 'string') {
            return {error: e}
        } else {
            return {error: 'Unknown error when fetching events over time '}
        }
    }
}

export const getDrugsFromEvents = async (searchTerm: string, searchOptions: SearchOptions) => {
    const url = generatePath([searchTerm], searchOptions, 'drug')
    try {
        const response = await fetch(url)
        const data = await response.json()
        return {result: processDrugEvents(data)}
    } catch (e: unknown) {
        if (e instanceof Error) {
            return {error: e.message}
        } else if (typeof e === 'string') {
            return {error: e}
        } else {
            return {error: 'Unknown error when searching by ADE'}
        }
    }
}

const demographicDataCache: {[key: string]: any} = {}

export const getSideEffectsForDemographics = async (drug: string) => {
    const demographicGroups: DemographicGroups[] = [
        { name: 'Young Boys', sex: '1', age: [0, 18], def: 'Males aged 0-18 years'},
        { name: 'Young Girls', sex: '2', age: [0, 18], def: 'Females aged 0-18 years'},
        { name: 'Men', sex: '1', age: [19, 59], def: 'Males aged 19-59 years'},
        { name: 'Women', sex: '2', age: [19, 59], def: 'Females aged 19-59 years'},
        { name: 'Elderly Men', sex: '1', age: [60, 120], def: 'Males aged 60-120 years'},
        { name: 'Elderly Women', sex: '2', age: [60, 120], def: 'Females aged 60-120 years'},
    ];

    // If the data for this drug is already in the cache, return it
    if (demographicDataCache[drug]) {
        return demographicDataCache[drug];
    }

    const results: {[key: string]: any} = {};

    for (const group of demographicGroups) {
        const searchOptions = {
            searchBy: {
                value: 'patient.drug.medicinalproduct',
                index: 9999,
                label: '',
                type: 'searchBy',
                enabled: true,
                param: 'medicinal_product'
            },
            sex: {
                value: `patient.patientsex:${group.sex}`,
                index: 0,
                label: '',
                type: 'sex',
                enabled: true
            },
            age: {
                min: {
                    value: group.age[0].toString(),
                    index: 0,
                    label: '',
                    type: 'age'
                },
                max: {
                    value: group.age[1].toString(),
                    index: 1,
                    label: '',
                    type: 'age'
                },
                enabled: true
            },
            country: {
                value: 'US',
                index: 0,
                label: 'United States',
                type: 'country',
                enabled: false
            }
        };

        const url = generatePath([drug], searchOptions, 'reaction')
        const response = await fetch(url)
        const data = await response.json()
        if(!response.ok) {
            break
        }
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
