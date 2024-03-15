import React from 'react';
import ReactApexChart from 'react-apexcharts';
import { ChartDataPoint, DemographicDataType, URLParams } from '../../types';
import useGeneralOptionsStore from '../../stores/generalOptionsStore';
import _ from 'lodash';
import { searchAgeGroups, searchSex, chartColors, oppositeAggregation } from '../../constants';
import useDemographicStore from '../../stores/demographicStore';
import { useUrlParams } from '../../hooks/useUrlParams';
import { useDemographicData as useFdaDemographicData } from '../../hooks/useDemographicData';
import { ApexOptions } from 'apexcharts';
import { capitalizeFirstLetter } from '../../utils/utils';
import { Row } from 'react-bootstrap';
import usePmDemographicData from '../../hooks/usePmDemographicData';

const getParamsArray = (term: string, searchBy: string, searchMode: string): URLParams[] => {
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

/* GROUP THE DATA BY SELECTED AGGREGATE TYPE (MALE, FEMALE, AGE GROUP ETC) */
const groupData = (data: DemographicDataType[], filterType: keyof DemographicDataType['params']) => {
  const aggregatedData: Record<string, DemographicDataType[]> = {};
  data.forEach((entry) => {
    const key = entry.params[filterType];
    if (key) {
      aggregatedData[key] = aggregatedData[key] || [];
      aggregatedData[key]!.push(entry);
    }
  });
  return aggregatedData;
};

/* ADVANCED VIEW DATA TRANSFORMATION */
const transformDataGranular = (data: DemographicDataType[], aggregateType: string) => {
  const aggregation = oppositeAggregation[aggregateType];

  if (!aggregation || !data) {
    return { series: [], labels: [] };
  }

  const labels = data.map((entry) => entry.params[aggregation]);
  const uniqueTerms = [...new Set(data.flatMap((item) => item.data.map((d) => d.x)))];

  const series = uniqueTerms.map((term) => ({
    name: capitalizeFirstLetter(term),
    data: data.map((item) => {
      const termData = item.data.find((d) => d.x === term);
      return termData ? termData.y : 0;
    }),
  }));

  return { series, labels };
};

/* SIMPLE VIEW DATA TRANSFORMATION */
const transformDataSimple = (data: DemographicDataType[], limit = 10) => {
  if (!data) {
    return { series: [], labels: [] };
  }

  // Combine all the data into one array of objects
  const combinedData = data.reduce((acc, curr) => acc.concat(curr.data), [] as ChartDataPoint[]);

  const combinedDataNoDupes = combinedData.reduce(
    (acc, { x, y }) => {
      acc[x] = (acc[x] || 0) + y;
      return acc;
    },
    {} as Record<string, number>
  );

  const combinedDataArray = Object.entries(combinedDataNoDupes)
    .map(([x, y]) => ({ x: capitalizeFirstLetter(x), y }))
    .sort((a, b) => b.y - a.y);

  return {
    series: [
      {
        name: 'Count',
        data: combinedDataArray
          .filter((entry) => entry.x !== 'Other')
          .slice(0, limit)
          .map((entry) => entry.y),
      },
    ],
    labels: combinedDataArray
      .filter((entry) => entry.x !== 'Other')
      .slice(0, limit)
      .map((entry) => entry.x),
  };
};

/* ALTERNATE BETWEEN FDA AND PM */
const useDynamicDemoData = (paramsArray: URLParams[], source: 'fda' | 'pm') => {
  return source === 'fda' ? useFdaDemographicData(paramsArray) : usePmDemographicData(paramsArray);
};

interface DemographicComparsionChartTypes {
  aggregateType: string;
  currentPageKey: string;
  onDataStatusChange: (status: boolean) => void;
  advancedView: boolean;
  source: 'fda' | 'pm';
}

const DemographicComparsionChart: React.FC<DemographicComparsionChartTypes> = ({
  aggregateType,
  currentPageKey,
  onDataStatusChange,
  advancedView,
  source,
}) => {
  const theme = useGeneralOptionsStore((state) => state.theme);

  const term = useDemographicStore((state) => state.demographicTerm);
  const searchBy = useDemographicStore((state) => state.demographicType);
  const searchMode = useUrlParams().params.searchMode;
  const [groupPageKeys, setGroupPageKeys] = useDemographicStore((state) => [
    state.groupPageKeys,
    state.setGroupPageKeys,
  ]);

  const paramsArray = getParamsArray(term, searchBy, searchMode);

  const { paramDataArray, isError } = useDynamicDemoData(paramsArray, source);

  const hasData = React.useMemo(
    () => !!(paramDataArray && !isError && paramDataArray.length !== 0),
    [paramDataArray, isError]
  );

  React.useEffect(() => {
    onDataStatusChange(hasData);
  }, [hasData]);

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

  if (!hasData) {
    return null;
  }

  const aggregatedDataForKey =
    source === 'fda' ? aggregatedData[currentPageKey] : aggregatedData[currentPageKey]?.slice(0, 10);

  const { series, labels } = advancedView
    ? transformDataGranular(aggregatedDataForKey || [], aggregateType)
    : transformDataSimple(aggregatedDataForKey || []); 

  const totalTermCount = series.reduce((acc, curr) => {
    return acc + curr.data.reduce((acc, curr) => acc + curr, 0);
  }, 0);

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
      stacked: advancedView ? true : false,
      stackType: '100%',
      background: theme === 'dark' ? '#212529' : '',
    },
    legend: {
      show: advancedView ? true : false,
    },
    plotOptions: {
      bar: {
        horizontal: true,
        distributed: advancedView ? false : true,
        barHeight: '60%',
        borderRadius: 8,
        borderRadiusWhenStacked: 'all',
        borderRadiusApplication: advancedView ? 'around' : 'end',
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
        style: {
          fontSize: '16px',
        },
      },
    },
    dataLabels: {
      enabled: true,
      formatter: (val: number) => {
        return advancedView ? `${val.toPrecision(3)}%` : `${((val / totalTermCount!) * 100).toPrecision(3)}%`;
      },
      style: {
        fontSize: '16px',
      },
    },
    tooltip: {
      x: {
        formatter: (val: number) => {
          return advancedView ? `Group: ${val}` : `${val}`;
        },
      },
      y: {
        formatter: (val: number) => {
          return val.toLocaleString();
        },
      },
    },
  };

  return (
    <div>
      <Row className={'text-center'}>
        <span className={'text-secondary'}>{totalTermCount?.toLocaleString()} terms for this chart in total</span>
      </Row>
      <ReactApexChart options={chartOptions} type="bar" series={series} />
    </div>
  );
};

export default DemographicComparsionChart;
