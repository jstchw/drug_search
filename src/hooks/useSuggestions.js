import {useEffect, useState} from "react";
import Fuse from "fuse.js";

export const useSuggestions = (apiEndpoint, fuseOptions) => {
    const [fuse, setFuse] = useState(null);

    useEffect(() => {
        async function fetchSuggestions() {
            try {
                const response = await fetch(apiEndpoint);
                const data = await response.json();
                setFuse(new Fuse(data, fuseOptions));
            } catch (error) {
                console.error('Error fetching suggestions:', error);
            }
        }

        fetchSuggestions();
    }, [apiEndpoint, fuseOptions]);

    return fuse;
}