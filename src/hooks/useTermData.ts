import {useUrlParams} from "./useUrlParams";
import {generatePath, processTermData} from "../utils/utils";
import {useEffect, useState} from "react";
import {ChartDataPoint, FDARawData, ResultItem} from "../types";

export const useTermData = () => {
    const { params } = useUrlParams();
    const url = generatePath(params);

    const [data, setData] = useState<ChartDataPoint[] | null>(null); // Set initial data to null
    const [error, setError] = useState<unknown | null>(null);
    const [loading, setLoading] = useState<boolean>(true);

    useEffect(() => {
        const fetchData = async () => {
            try {
                const response = await fetch(url);
                // Check if the response is ok (status in the range 200-299)
                if (!response.ok) {
                    setError(`HTTP error: ${response.status}`);
                }
                const result: FDARawData = await response.json();
                const processedData = processTermData(result.results as ResultItem[]);
                setData(processedData); // Update the state with the parsed result
                setLoading(false); // Set loading to false after data is fetched
                setError(null);
            } catch (e: unknown) {
                setError(e); // Set error if an exception occurs
                setLoading(false); // Ensure loading is set to false on error
            }
        };

        fetchData();
    }, [url]); // Dependency array with url to re-run effect when url changes

    return {data, error, loading}
}