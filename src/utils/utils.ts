import {ChartDataPoint, ResultItem, SearchTypesMap, TimeEventData, URLParams, SearchOptionsType} from "../types";
import {baseFdaUrl, searchTypes} from "../constants";


export const processTermData = (data: ResultItem[]): ChartDataPoint[] => {
    return data.map(item => ({
        x: item.term,
        y: item.count,
    }))
}

export const processYearData = (data: TimeEventData[]): ChartDataPoint[] => {
    const yearlyData = data.reduce((acc: { [year: string]: number }, entry: TimeEventData) => {
        const year = entry.time.substring(0, 4) // Accessing the year from the time string (e.g. 2019-01-01)
        acc[year] = (acc[year] || 0) + entry.count // If the year exists in the accumulator, add the count to it, otherwise set it to 0 and add the count to it (this is to avoid undefined errors)
        return acc
    }, {})

    const sortedEntries: [string, number][] = Object.entries(yearlyData).sort((a, b) => a[0].localeCompare(b[0]))

    return sortedEntries.map(([time, count]): ChartDataPoint => ({
        x: time,
        y: count,
    }))
}

const formatDate = (date: Date) => {
    const year = date.getFullYear() - 1
    const month = date.getMonth() + 1
    const day = date.getDate()
    return `${year}${month < 10 ? `0${month}` : month}${day < 10 ? `0${day}` : day}`
}

export const generatePath = (params: URLParams, countType?: string) => {
    const fromDate = `20040101`
    const limit = 50
    const toDate = formatDate(new Date())
    const whatToCount: SearchTypesMap = {
        generic_name: 'patient.reaction.reactionmeddrapt.exact',
        brand_name: 'patient.reaction.reactionmeddrapt.exact',
        receiveDate: 'receivedate',
        side_effect: 'patient.drug.activesubstance.activesubstancename.exact'
    }

    const searchParts = [`(receivedate:[${fromDate}+TO+${toDate}])`]

    // After the Search page has been translated to TS, work on the searchBy in params
    // It needs to be passed as a whole object, not just the value and then the value needs to be extracted
    if (params.terms) {
        const encodedTerms = params.terms.flat().join('+AND+')
        searchParts.push(`(${mapParamToValue(params.searchBy, searchTypes)}:${encodedTerms})`)
    }

    if (params.sex) {
        searchParts.push(params.sex)
    }

    if (params.age && params.age.min && params.age.max) {
        searchParts.push(`patient.patientonsetage:[${params.age.min}+TO+${params.age.max}]`)
    }

    if (params.country) {
        searchParts.push(`occurcountry:"${params.country}"`)
    }

    if (countType) {
        return `${baseFdaUrl}?search=${searchParts.join('+AND+')}&count=${whatToCount[countType]}`
    } else {
        return `${baseFdaUrl}?search=${searchParts.join('+AND+')}&count=${whatToCount[params.searchBy]}&limit=${limit}`
    }
}

export const mapParamToLabel = (param: string, options: SearchOptionsType[]): string => {
    const option = options.find(option => option.param === param)
    return option ? option.label : 'Not Specified'
}

export const mapParamToValue = (param: string, options: SearchOptionsType[]): string => {
    const option = options.find(option => option.param === param)
    return option ? option.value : ''
}