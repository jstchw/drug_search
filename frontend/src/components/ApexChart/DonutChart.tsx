import Chart from 'react-apexcharts';
import useGeneralOptionsStore from '../../stores/generalOptionsStore';
import React from 'react';
import { generatePath, fetchData, processTermData } from '../../utils/utils';
import useDemographicStore from '../../stores/demographicStore';
import { useUrlParams } from '../../hooks/useUrlParams';
import { useQuery } from 'react-query';
import { BackendDataType, ChartDataPoint, FDARawData, ResultItem } from '../../types';
import { ageGroupsFDA, sexGroupsFDA, backendUrl, chartColors } from '../../constants';
import _ from 'lodash';
import { ApexOptions } from 'apexcharts';

type RadarChartReturnType = {
  data: ChartDataPoint[] | undefined;
  isError: boolean;
  isLoading: boolean;
};

const fetchFdaDistributionData = (type: string): RadarChartReturnType => {
  const term = useDemographicStore((state) => state.demographicTerm);
  const searchBy = useDemographicStore((state) => state.demographicType);
  const {
    params: { searchMode },
  } = useUrlParams();
  // Mimicking params
  const params = { terms: [term], searchBy, searchMode };
  const url = generatePath(params, type, 10);

  const { data, isError, isLoading } = useQuery(
    ['distributionChart', url],
    () => fetchData(url) as Promise<FDARawData>,
    {
      staleTime: 3600000,
      retry: false,
      select: (data) => processTermData(data.results as ResultItem[]),
    }
  );

  return { data, isError, isLoading };
};

const fetchPmdDistributionData = (type: string): RadarChartReturnType => {
  const term = useDemographicStore((state) => state.demographicTerm);
  const searchBy = useDemographicStore((state) => state.demographicType);
  const {
    params: { searchMode },
  } = useUrlParams();

  let url = '';

  if (type === 'age_group') {
    url =
      `${backendUrl}/drug/get_pm_age_distribution?` +
      `terms=${term}&` +
      `search_mode=${searchMode}&` +
      `search_type=${searchBy}&`;
  } else if (type === 'patient_sex') {
    url =
      `${backendUrl}/drug/get_pm_sex_distribution?` +
      `terms=${term}&` +
      `search_mode=${searchMode}&` +
      `search_type=${searchBy}&`;
  }

  const { data, isError, isLoading } = useQuery(
    ['distributionChart', url],
    () => fetchData(url) as Promise<BackendDataType>,
    {
      staleTime: 3600000,
      retry: false,
      select: (data) => data.data as ChartDataPoint[],
    }
  );

  return { data, isError, isLoading };
};

const useDynamicDistributionData = (source: string, type: string) => {
  return source === 'fda' ? fetchFdaDistributionData(type) : fetchPmdDistributionData(type);
};

const augmentDistributionData = (data: ChartDataPoint[], type: string) => {
  const dataLabels = data.map((entry) => {
    switch (type) {
      case 'age_group':
        return ageGroupsFDA[entry.x] || entry.x;
      case 'patient_sex':
        return sexGroupsFDA[entry.x] || entry.x;
      default:
        return entry.x;
    }
  });

  const dataSeries = data.map((entry) => entry.y);

  return { dataLabels, dataSeries };
};

interface DonutChartProps {
  source: string;
  type: string;
  onDataStatusChange: (status: boolean) => void;
}

const DonutChart: React.FC<DonutChartProps> = ({ source, type, onDataStatusChange }) => {
  const theme = useGeneralOptionsStore((state) => state.theme);

  const { data, isError } = useDynamicDistributionData(source, type);

  const hasData = React.useMemo(() => !!(data && !isError && data.length !== 0), [data, isError]);

  React.useEffect(() => {
    onDataStatusChange(hasData);
  }, [hasData]);

  if (!data || isError) {
    return;
  }

  const { dataLabels, dataSeries } = augmentDistributionData(data, type);

  const reportCount = dataSeries.reduce((acc, curr) => acc + curr, 0);

  const chartOptions: ApexOptions = {
    theme: {
      mode: theme,
    },
    labels: dataLabels,
    legend: {
      show: true,
      position: 'bottom',
    },
    tooltip: {
      y: {
        formatter: (val: number) => {
          return val.toLocaleString();
        },
      },
      style: {
        fontSize: '16px',
      },
    },
    colors: chartColors,
    chart: {
      toolbar: {
        show: false,
        tools: {
          zoom: false,
          zoomin: false,
          zoomout: false,
        },
      },
      background: theme === 'dark' ? '#212529' : '',
    },
  };

  return (
    <div className={'d-flex flex-column'}>
      {type === 'age_group' ? (
        <span className={'fs-4 fw-light mb-1'}>Age distribution</span>
      ) : (
        <span className={'fs-4 fw-light mb-1'}>Sex distribution</span>
      )}
      <span className={'text-secondary mb-2'}>From {reportCount.toLocaleString()} reports</span>
      <Chart options={chartOptions} series={dataSeries} type="donut" />
    </div>
  );
};

export default DonutChart;
