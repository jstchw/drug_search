import Chart from "react-apexcharts";
import { getDonutChartOptions } from "./chartOptions";
import { ThemeContext } from "../../contexts/ThemeContext";
import React from "react";
import { generatePath, fetchData, processTermData } from "../../utils/utils";
import useDemographicStore from "../../stores/demographicStore";
import { useUrlParams } from "../../hooks/useUrlParams";
import { useQuery } from "react-query";
import { ChartDataPoint, ResultItem } from "../../types";
import { ageGroupsFDA, sexGroupsFDA } from "../../constants";
import _ from "lodash";

type RadarChartReturnType = {
  data: ChartDataPoint[] | undefined;
  isError: boolean;
  isLoading: boolean;
};

const fetchTreeMapData = (type: string): RadarChartReturnType => {
  const term = useDemographicStore((state) => state.demographicTerm);
  const searchBy = useDemographicStore((state) => state.demographicType);
  const {
    params: { searchMode },
  } = useUrlParams();
  // Mimicking params
  const params = { terms: [term], searchBy, searchMode };
  const url = generatePath(params, type, 10);

  const { data, isError, isLoading } = useQuery(
    ["treeMapChart", url],
    () => fetchData(url),
    {
      staleTime: 3600000,
      retry: false,
      select: (data) => processTermData(data.results as ResultItem[]),
    },
  );

  return { data, isError, isLoading };
};

const augmentDonutData = (data: ChartDataPoint[], type: string) => {
  const dataLabels = data.map((entry) => {
    switch (type) {
      case "age_group":
        return ageGroupsFDA[entry.x] || entry.x;
      case "patient_sex":
        return sexGroupsFDA[entry.x] || entry.x;
      default:
        return entry.x;
    }
  });

  const dataSeries = data.map((entry) => entry.y);

  return { dataLabels, dataSeries };
};

const DonutChart = ({ type }: { type: string }) => {
  const { theme } = React.useContext(ThemeContext);
  const { data, isError } = fetchTreeMapData(type);

  if (!data || isError) {
    return;
  }

  const { dataLabels, dataSeries } = augmentDonutData(data, type);

  const reportCount = dataSeries.reduce((acc, curr) => acc + curr, 0);

  const chartOptions = getDonutChartOptions(theme);

  const specificOptions = {
    labels: dataLabels,
  };

  const options = _.merge(chartOptions, specificOptions);

  return (
    <div className={"d-flex flex-column"}>
      <span className={"text-secondary mb-2"}>
        From {reportCount.toLocaleString()} reports
      </span>
      <Chart options={options} series={dataSeries} type="donut" />
    </div>
  );
};

export default DonutChart;
