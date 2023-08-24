import Cookies from "js-cookie";
import React from "react";
import { SearchOptions } from "../types";

const updateSearchHistory = (
    capitalizedTerms: string[],
    searchHistory: string[][] | null,
    setSearchHistory: React.Dispatch<React.SetStateAction<string[][] | null>>,
    searchOptions: SearchOptions,
) => {
    console.log('searchOptions', searchOptions);

    searchHistory && searchHistory[0] !== undefined && console.log(compareDiffOrderStrings(searchHistory[0], capitalizedTerms));

    if (searchHistory && searchHistory[0] !== undefined) {
        // Finding the index of the duplicate search term
        const duplicateIndex = searchHistory.findIndex(
            (search: string[]) => compareDiffOrderStrings(search, capitalizedTerms)
        );

        // If a duplicate is found, remove it
        if (duplicateIndex !== -1) {
            searchHistory.splice(duplicateIndex, 1);
        }

        // Add the ORIGINAL capitalizedTerms and update the state
        setSearchHistory([capitalizedTerms, ...searchHistory.slice(0, 4)]);
        Cookies.set('searchHistory', JSON.stringify([capitalizedTerms, ...searchHistory.slice(0, 4)]));
    } else {
        setSearchHistory([capitalizedTerms]);
        Cookies.set('searchHistory', JSON.stringify([capitalizedTerms]));
    }
}

export default updateSearchHistory;

const compareDiffOrderStrings = (str1: string[], str2: string[]) => {
    const sortedStr1 = str1.slice().sort().join();
    const sortedStr2 = str2.slice().sort().join();

    return sortedStr1 === sortedStr2;
}
