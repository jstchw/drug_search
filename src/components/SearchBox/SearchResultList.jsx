import React from "react";
import Chart from "react-apexcharts";

const SearchResultList = ({ searchResults }) => {
    return (
        <div>
            {searchResults.results && searchResults.results.map((result, index) => (
                <div key={index}>
                    <p>Term: {result.term}</p>
                    <p>Count: {result.count}</p>

                </div>
            ))}
        </div>
    );
}

export default SearchResultList;