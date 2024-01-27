import { useUrlParams } from "./useUrlParams";
import { generatePath, processTermData } from "../utils/utils";
import { useEffect, useState } from "react";
import { ChartDataPoint, FDARawData, ResultItem } from "../types";

export const useTermData = () => {
  const { params } = useUrlParams();
  const url = generatePath(params);

  const [data, setData] = useState<ChartDataPoint[] | null>(null);
  const [error, setError] = useState<unknown | boolean>(false);
  const [loading, setLoading] = useState<boolean>(true);

  useEffect(() => {
    setError(false);
    let isMounted = true; // Flag to check if component is still mounted

    const fetchData = async () => {
      localStorage.setItem("isFetching", "true"); // Set fetching flag
      try {
        const response = await fetch(url);
        if (!response.ok) {
          setError(`HTTP error: ${response.status}`);
        }
        const result: FDARawData = await response.json();
        const processedData = processTermData(result.results as ResultItem[]);
        if (isMounted) {
          setData(processedData);
          setError(null);
        }
        localStorage.setItem(
          "termData",
          JSON.stringify({ processedData, params, timestamp: Date.now() }),
        );
      } catch (e: unknown) {
        if (isMounted) {
          setError(e);
        }
      } finally {
        if (isMounted) {
          setLoading(false);
        }
        localStorage.removeItem("isFetching"); // Clear fetching flag
      }
    };

    const waitForOtherFetch = async () => {
      while (localStorage.getItem("isFetching")) {
        await new Promise((resolve) => setTimeout(resolve, 100));
      }
      if (isMounted) {
        fetchData();
      }
    };

    const rawCachedData = localStorage.getItem("termData");
    if (rawCachedData) {
      const {
        processedData: cachedData,
        params: cachedParams,
        timestamp,
      } = JSON.parse(rawCachedData);
      if (
        Date.now() - timestamp < 1000 * 60 * 60 * 24 &&
        JSON.stringify(cachedParams) === JSON.stringify(params)
      ) {
        setData(cachedData);
        setLoading(false);
      } else if (localStorage.getItem("isFetching")) {
        waitForOtherFetch();
      } else {
        fetchData();
      }
    } else if (localStorage.getItem("isFetching")) {
      waitForOtherFetch();
    } else {
      fetchData();
    }

    return () => {
      isMounted = false; // Clean up to avoid setting state on unmounted component
    };
  }, [params, url]);

  return { data, error, loading };
};
