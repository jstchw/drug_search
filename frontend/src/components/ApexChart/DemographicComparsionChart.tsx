import React from 'react';
import ReactApexChart from 'react-apexcharts';
import useGeneralOptionsStore from '../../stores/generalOptionsStore';
import _ from 'lodash';
import { searchAgeGroups, searchSex, chartColors } from '../../constants';
import useDemographicStore from '../../stores/demographicStore';
import { useUrlParams } from '../../hooks/useUrlParams';
import { ApexOptions } from 'apexcharts';
import { Row } from 'react-bootstrap';
import useDemographicData from '../../hooks/useDemographicData';
import { BarLoader } from 'react-spinners';
import { motion } from 'framer-motion';
import { generateRequestArgs, valueToPercentage } from '../../utils/utils';

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
  const searchType = useDemographicStore((state) => state.demographicType);
  const searchMode = useUrlParams().params.searchMode;
  const [groupPageKeys, setGroupPageKeys] = useDemographicStore((state) => [
    state.groupPageKeys,
    state.setGroupPageKeys,
  ]);


  const requestArgs = generateRequestArgs(term, searchType, aggregateType, currentPageKey, advancedView, searchMode);

  const { data, isError, isLoading } = useDemographicData(requestArgs, source);

  const hasData = React.useMemo(
    () => !!(data && !isError),
    [data, isError]
  );

  React.useEffect(() => {
    onDataStatusChange(true);
  }, [hasData]);

  React.useEffect(() => {
    let newKeys;
    if (aggregateType === 'Sex') {
      newKeys = Object.values(searchSex).map((sexOption) => sexOption.label);
    } else if (aggregateType === 'Age') {
      newKeys = Object.keys(searchAgeGroups);
    }

    if (!_.isEqual(groupPageKeys, newKeys ?? [])) {
      setGroupPageKeys(newKeys ?? []);
    }
  }, [data, groupPageKeys, setGroupPageKeys]);

  if (!hasData) {
    return null;
  }

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
      categories: data?.categories,
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
      formatter: function(val, opts) {
        // Assuming your raw data is accessible, e.g., in a `series` array
        const rawValue = data.series[opts.seriesIndex].data[opts.dataPointIndex];
        return rawValue?.toLocaleString(); // This will display raw values instead of percentages
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
        formatter: (val: number, opts) => {
          if (advancedView) {
            const categoryIndex = opts.dataPointIndex;
            const sumForCategory = data.series.reduce((acc, series) => acc + series.data[categoryIndex], 0);
            return `${valueToPercentage(val, sumForCategory).toPrecision(3)}% out of ${sumForCategory.toLocaleString()}`;
          } else {
            return `${valueToPercentage(val, data?.total_count).toPrecision(3)}% out of ${data?.total_count.toLocaleString()}`;
          }
        },
      },
    },
  };

  return (
    <>
    {isLoading && (
      <div className={'d-flex justify-content-center mb-4'}>
        <BarLoader className={'text-red'} loading={isLoading} speedMultiplier={2}/>
      </div>
    )}
    <motion.div layout transition={{ type: 'spring', stiffness: 260, damping: 20, delay: 0.1 }}
      style={{
        filter: isLoading ? 'blur(0.5em) grayscale(1)' : 'none',
        pointerEvents: isLoading ? 'none' : 'auto',
        }}
    >
      {data?.total_count && (
      <Row className={'text-center'}>
        <span className={'text-secondary'}>Percentage is calculated using only the presented data</span>
        <span className={'text-secondary'}>Some data may not be displayed due to low significance</span>
        <span className={'text-secondary mt-2'}>{data.total_count.toLocaleString()} events collected for {data.categories.length} {advancedView ? 'categories' : 'terms'}</span>
      </Row>
      )}
      <ReactApexChart options={chartOptions} series={data?.series} type="bar" />
    </motion.div>
    </>
  );
};

export default DemographicComparsionChart;
