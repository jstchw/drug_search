import React from 'react';
import { useDynamicTimeSeries } from '../../hooks/useDynamicTimeSeries';
import { ApexOptions } from 'apexcharts';
import ReactApexChart from 'react-apexcharts';
import { motion } from 'framer-motion';
import useGeneralOptionsStore from '../../stores/generalOptionsStore';
import { Row } from 'react-bootstrap';

interface CasesTimeSeriesChartProps {
  noFilterRequest?: boolean;
  onRender: () => void;
  isPubmed?: boolean;
}

const CasesTimeSeriesChart: React.FC<CasesTimeSeriesChartProps> = ({
  noFilterRequest = false,
  onRender,
  isPubmed = false,
}) => {
  const theme = useGeneralOptionsStore((state) => state.theme);

  const { timeSeriesData, timeSeriesCount, isError } = useDynamicTimeSeries(isPubmed, noFilterRequest);

  React.useEffect(() => {
    if (timeSeriesData && !isError) {
      onRender();
    }
  }, [timeSeriesData, isError]);

  if (!timeSeriesData || isError) {
    return null;
  }

  const chartData = {
    series: [
      {
        name: 'Reports',
        data: timeSeriesData,
      },
    ],
  };

  const options: ApexOptions = {
    theme: {
      mode: theme,
    },
    colors: ['#59768A'],
    legend: {
      show: true,
      position: 'bottom',
    },
    chart: {
      type: 'area',
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
    stroke: {
      curve: 'smooth',
    },
    fill: {
      type: 'gradient',
    },
    dataLabels: {
      enabled: false,
    },
    xaxis: {
      type: 'category',
      labels: {
        style: {
          fontSize: '14px',
          colors: theme === 'dark' ? '#ACB5BD' : '',
        },
      },
      tooltip: {
        enabled: false,
      },
    },
    yaxis: {
      labels: {
        style: {
          fontSize: '14px',
          colors: theme === 'dark' ? '#ACB5BD' : '',
        },
        formatter: (value: number) => {
          return value.toLocaleString();
        },
      },
    },
  };

  return (
    <>
      {timeSeriesCount !== undefined && (
        <Row className={'text-center'}>
          <span className={'text-secondary'}>{timeSeriesCount.toLocaleString()} reports in total</span>
        </Row>
      )}
      <motion.div
        layout
        initial={{ opacity: 0 }}
        animate={{ opacity: 1 }}
        transition={{
          type: 'spring',
          stiffness: 260,
          damping: 20,
        }}
        exit={{ opacity: 0 }}
      >
        <ReactApexChart options={options} series={chartData.series} type={options.chart?.type} />
      </motion.div>
    </>
  );
};

export default CasesTimeSeriesChart;
