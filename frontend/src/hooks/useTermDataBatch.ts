import { FDARawData, URLParams } from "../types";
import { useEffect, useState } from "react";
import { ChartDataPoint, ResultItem } from "../types";
import { generatePath, processTermData } from "../utils/utils";

export const useTermDataBatch = (paramsArray: URLParams[]) => {
  const [paramDataArray, setParamDataArray] = useState<
    { params: URLParams; data: ChartDataPoint[] }[] | null
  >(null);
  const [error, setError] = useState<boolean>(false);
  const [loading, setLoading] = useState<boolean>(true);

  useEffect(() => {
    setError(false);
    let isMounted = true; // Flag to check if component is still mounted

    const fetchData = async () => {
      try {
        const fetchPromises = paramsArray.map((params) => {
          const url = generatePath(params, undefined, 10);
          return fetch(url)
            .then((response) => {
              if (!response.ok)
                throw new Error(`HTTP error: ${response.status}`);
              return response.json();
            })
            .then((result: FDARawData) => ({
              params,
              data: processTermData(result.results as ResultItem[]),
            }));
        });

        const results = await Promise.all(fetchPromises);

        if (isMounted) {
          setParamDataArray(results);
          setError(false);
        }
      } catch (e: unknown) {
        if (isMounted) {
          setError(true);
        }
      } finally {
        if (isMounted) {
          setLoading(false);
        }
      }
    };

    void fetchData();

    return () => {
      isMounted = false;
    };
  }, []);

  return { paramDataArray, error, loading };
};
