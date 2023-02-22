import React from "react";

const SearchResultList = ({ searchResults }) => {
    return (
        <div>
            {searchResults?.results?.map((result) => (
                <div key={result.safetyreportid}>
                    <p>Receive Date: {result.receivedate}</p>
                </div>
            ))}
        </div>
    );
}

export default SearchResultList;