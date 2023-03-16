import React from "react";

const SearchResultList = ({ searchResults }) => {
    return (
        <div>
            {/*{searchResults?.results?.map((result) => (*/}
            {/*    <div key={result.safetyreportid}>*/}
            {/*        <p>Receive Date: {result.receivedate}</p>*/}
            {/*    </div>*/}
            {/*))}*/}
            <pre>
                {JSON.stringify(searchResults, null, 2)}
            </pre>
        </div>
    );
}

export default SearchResultList;