import { useEffect, useState } from "react";
import { generatePath } from "../utils/utils";
import { useUrlParams } from "./useUrlParams";
import {ChartDataPoint, FDARawData, TimeEventData} from "../types";
import {processYearData} from "../utils/utils";

export const useTimeSeriesData = () => {
    const { params } = useUrlParams();
    const url = generatePath(params, 'receiveDate');

    const [data, setData] = useState<ChartDataPoint[] | null>(null);
    const [loading, setLoading] = useState<boolean>(true);
    const [error, setError] = useState<unknown | boolean>(false);

    useEffect(() => {
        setError(false); // Reset error
        const fetchData = async () => {
            try {
                const response = await fetch(url);
                // Check if the response is ok (status in the range 200-299)
                if (!response.ok) {
                    setError(`HTTP error: ${response.status}`);
                }
                const result: FDARawData = await response.json();
                const processedData = processYearData(result.results as TimeEventData[]);
                setData(processedData); // Update the state with the parsed result
                setLoading(false); // Set loading to false after data is fetched
            } catch (e: unknown) {
                setError(e); // Set error if an exception occurs
                setLoading(false); // Ensure loading is set to false on error
            }
        };

        fetchData();
    }, [url]); // Dependency array with url to re-run effect when url changes

    return { data, loading, error };
};
