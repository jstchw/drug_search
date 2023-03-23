import React from "react";
import PieChart from "../PieChart/PieChart";

const SearchResultObject = ({ searchResults }) => {
    const { termCountDict, totalCount } = searchResults.result

    return (
        <div>
            {searchResults && <PieChart termCountDict={termCountDict} totalCount={totalCount} />}
        </div>
    );
}

export default SearchResultObject;