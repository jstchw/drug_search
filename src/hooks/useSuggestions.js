import {useEffect, useState} from "react";
import Fuse from "fuse.js";

export const useSuggestions = (apiEndpoint, fuseOptions, searchBy) => {
    const [fuse, setFuse] = useState(null);

    useEffect(() => {
        async function fetchSuggestions() {
            try {
                // Form a URL based on the searchBy parameter and fetch the data
                const url = `${apiEndpoint}?searchBy=${encodeURIComponent(searchBy)}`;
                const response = await fetch(url);
                const data = await response.json();
                setFuse(new Fuse(data, fuseOptions));
            } catch (error) {
                console.error('Error fetching suggestions:', error);
            }
        }

        fetchSuggestions();
    }, [apiEndpoint, fuseOptions, searchBy]);

    return fuse;
}