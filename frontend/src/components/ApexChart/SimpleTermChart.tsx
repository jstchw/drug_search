import React from "react";
import ReactApexChart from "react-apexcharts";
import { DemographicData, URLParams } from "../../types";
import { ThemeContext } from "../../contexts/ThemeContext";
import _ from "lodash";
import {
  searchAgeGroups,
  searchSex,
  chartColors,
  oppositeAggregation,
} from "../../constants";
import useDemographicStore from "../../stores/demographicStore";
import { useUrlParams } from "../../hooks/useUrlParams";
import { useDemographicData } from "../../hooks/useDemographicData";
import { ApexOptions } from "apexcharts";
import { capitalizeFirstLetter } from "../../utils/utils";

const getParamsArray = (
  term: string,
  searchBy: string,
  searchMode: string,
): URLParams[] => {
  const paramList: URLParams[] = [];

  Object.entries(searchAgeGroups).forEach(([_, ageGroup]) => {
    searchSex.forEach((sexOption) => {
      const params: URLParams = {
        terms: [term],
        searchBy,
        searchMode,
        sex: sexOption.param,
        age: {
          min: ageGroup.min.toString(),
          max: ageGroup.max.toString(),
        },
      };
      paramList.push(params);
    });
  });

  return paramList;
};

const groupData = (
  data: DemographicData[],
  filterType: keyof DemographicData["params"],
) => {
  const aggregatedData: Record<string, DemographicData[]> = {};

  data.forEach((entry) => {
    const key = entry.params[filterType];
    if (key) {
      aggregatedData[key] = aggregatedData[key] || [];
      aggregatedData[key]!.push(entry);
    }
  });

  return aggregatedData;
};

const transformData = (data: DemographicData[], aggregateType: string) => {
  const aggregation = oppositeAggregation[aggregateType];

  if (!aggregation || !data) {
    return { series: [], labels: [] };
  }

  const labels = data.map((entry) => entry.params[aggregation]);
  const uniqueTerms = [
    ...new Set(data.flatMap((item) => item.data.map((d) => d.x))),
  ];

  const series = uniqueTerms.map((term) => ({
    name: capitalizeFirstLetter(term),
    data: data.map((item) => {
      const termData = item.data.find((d) => d.x === term);
      return termData ? termData.y : 0;
    }),
  }));

  return { series, labels };
};

interface SimpleTermChartProps {
  aggregateType: string;
  currentPageKey: string;
}

const SimpleTermChart: React.FC<SimpleTermChartProps> = ({
  aggregateType,
  currentPageKey,
}) => {
  const { theme } = React.useContext(ThemeContext);

  const term = useDemographicStore((state) => state.demographicTerm);
  const searchBy = useDemographicStore((state) => state.demographicType);
  const searchMode = useUrlParams().params.searchMode;
  const [groupPageKeys, setGroupPageKeys] = useDemographicStore((state) => [
    state.groupPageKeys,
    state.setGroupPageKeys,
  ]);

  const paramsArray = getParamsArray(term, searchBy, searchMode);

  const { paramDataArray, isError } = useDemographicData(paramsArray);

  const aggregatedData = React.useMemo(() => {
    if (!paramDataArray || isError) {
      return {};
    }
    return groupData(paramDataArray, aggregateType);
  }, [paramDataArray, isError, aggregateType]);

  React.useEffect(() => {
    const newKeys = Object.keys(aggregatedData);

    if (!_.isEqual(groupPageKeys, newKeys)) {
      setGroupPageKeys(newKeys);
    }
  }, [aggregatedData, groupPageKeys, setGroupPageKeys]);

  const aggregatedDataForKey = aggregatedData[currentPageKey];
  const { series, labels } = transformData(
    aggregatedDataForKey || [],
    aggregateType,
  );

  const chartOptions: ApexOptions = {
    colors: chartColors,
    theme: {
      mode: theme,
    },
    chart: {
      toolbar: {
        show: false,
        tools: {
          zoom: false,
          zoomin: false,
          zoomout: false,
        },
      },
      stacked: true,
      stackType: "100%",
      background: theme === "dark" ? "#212529" : "",
    },
    legend: {
      show: true,
    },
    plotOptions: {
      bar: {
        horizontal: true,
        distributed: false,
        barHeight: "50%",
        borderRadius: 0,
        borderRadiusWhenStacked: "last",
      },
    },
    xaxis: {
      categories: labels,
      labels: {
        show: false,
      },
    },
    yaxis: {
      labels: {
        show: true,
      },
    },
    tooltip: {
      x: {
        formatter: (val: number) => {
          return `Group: ${val.toLocaleString()}`;
        },
      },
    },
  };

  return (
    <>
      <div className={"mb-3"}>
        <ReactApexChart options={chartOptions} type="bar" series={series} />
      </div>
    </>
  );
};

export default SimpleTermChart;
