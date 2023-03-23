import React from "react";
import PieChart from "../PieChart/PieChart";

const SearchResultObject = ({ searchResults, searchTerm }) => {
    const { termCountDict, totalCount } = searchResults.result

    return (
        <div>
            <h3>Search Results for {searchTerm}</h3>
            {searchResults && <PieChart termCountDict={termCountDict} totalCount={totalCount} />}
        </div>
    );
}

export default SearchResultObject;