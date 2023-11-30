import Cookies from "js-cookie";
import React from "react";
import { SearchOptions } from "../types";
import { SearchHistoryElement } from "../types";

const updateSearchHistory = (
    capitalizedTerms: string[],
    searchHistory: SearchHistoryElement[] | null,
    setSearchHistory: React.Dispatch<React.SetStateAction<SearchHistoryElement[] | null>>,
    searchOptions: SearchOptions,
) => {
    if (searchHistory && searchHistory[0] !== undefined) {
        // Finding the index of the duplicate search term
        const duplicateIndex = searchHistory.findIndex(
            (search: SearchHistoryElement) => compareDiffOrderStrings(search, capitalizedTerms)
        );

        // If a duplicate is found, remove it
        if (duplicateIndex !== -1) {
            searchHistory.splice(duplicateIndex, 1);
        }

        const newElement = { terms: capitalizedTerms, options: searchOptions };
        const newSearchHistory = [newElement, ...searchHistory.slice(0, 4)];

        // Add the ORIGINAL capitalizedTerms and update the state
        setSearchHistory(newSearchHistory);
        Cookies.set('searchHistory', JSON.stringify(newSearchHistory));
    } else {
        const newElement = { terms: capitalizedTerms, options: searchOptions };
        setSearchHistory([newElement]);
        Cookies.set('searchHistory', JSON.stringify([newElement]));
    }
}

export default updateSearchHistory;

const compareDiffOrderStrings = (item1: SearchHistoryElement, item2: string[]) => {
    const sortedStr1 = item1.terms.slice().sort().join();
    const sortedStr2 = item2.slice().sort().join();

    return sortedStr1 === sortedStr2;
}
