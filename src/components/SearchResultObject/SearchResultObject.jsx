import React from "react";
import PieChart from "../PieChart/PieChart";

const SearchResultObject = ({ searchResults, searchTerm }) => {
    const { termCountDict, totalCount } = searchResults.result

    return (
        <div>
            <h2>Search Results for {searchTerm}</h2>
            {searchResults && <PieChart termCountDict={termCountDict} totalCount={totalCount} />}
        </div>
    );
}

export default SearchResultObject;