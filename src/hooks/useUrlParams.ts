import {useLocation} from 'react-router-dom';
import React, {useMemo, useState} from 'react';
import {searchSex, searchTypes} from "../constants";
import {URLParams} from "../types";

type ValidationResult = {
    terms: boolean;
    searchBy: boolean;
    sex: boolean | null;
    min_age: boolean | null;
    max_age: boolean | null;
    country: boolean | null;
}

const isValidTerm = (term: string | null): boolean => !!term && term.length > 0;
const isValidSearchType = (type: string | null): boolean => !!type && searchTypes.some(t => t.param === type);
const isValidSex = (sex: string | null): boolean => !!sex && searchSex.some(s => s.param === sex);
const isValidAge = (age: string | null): boolean => age !== null && !isNaN(Number(age)) && Number(age) >= 0;

const isValidCountry = (country: string | null): boolean => country !== null;

const validateParams = (params: URLSearchParams): ValidationResult => {
    return {
        terms: isValidTerm(params.get('query')),
        searchBy: isValidSearchType(params.get('searchBy')),
        sex: params.has('sex') ? isValidSex(params.get('sex')) : null,
        min_age: params.has('min_age') ? isValidAge(params.get('min_age')) : null,
        max_age: params.has('max_age') ? isValidAge(params.get('max_age')): null,
        country: params.has('country') ? isValidCountry(params.get('country')) : null
    }
}

const parseUrlParams = (paramString: string, setError: React.Dispatch<boolean>) => {
    const params = new URLSearchParams(paramString)
    const validation = validateParams(params)

    const areParamsValid = Object.values(validation).every(val => val !== false)
    if(!areParamsValid) {
        setError(true)
        return {} as URLParams
    }
    const query = params.get('query')
    const terms = query ? query.split('-') : []
    return {
        terms: terms,
        searchBy: params.get('searchBy') as string,
        sex: params.get('sex'),
        age: {
            min: params.get('min_age'),
            max: params.get('max_age')
        },
        country: params.get('country')
    } as URLParams
}

export const useUrlParams = () => {
    const location = useLocation()
    const [error, setError] = useState<boolean>(false)



    const params = useMemo(() => parseUrlParams(location.search, setError), [location.search])

    return { params, error }
}