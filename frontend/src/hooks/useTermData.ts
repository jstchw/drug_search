import { useUrlParams } from "./useUrlParams";
import { generatePath, processTermData } from "../utils/utils";
import { useEffect, useState } from "react";
import { ChartDataPoint, FDARawData, ResultItem } from "../types";

export const useTermData = (noFilterRequest = false) => {
  let { params } = useUrlParams();

  if (noFilterRequest) {
    params = {
      ...params,
      sex: null,
      age: {
        min: null,
        max: null,
      },
      country: null,
    };
  }

  const url = generatePath(params);

  const [data, setData] = useState<ChartDataPoint[] | null>(null);
  const [error, setError] = useState<unknown | boolean>(false);
  const [loading, setLoading] = useState<boolean>(true);

  useEffect(() => {
    setError(false);
    let isMounted = true; // Flag to check if component is still mounted

    const fetchData = async () => {
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
      } else {
        void fetchData();
      }
    } else {
      void fetchData();
    }

    return () => {
      isMounted = false; // Clean up to avoid setting state on unmounted component
    };
  }, [url]);

  return { data, error, loading };
};
