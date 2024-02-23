import React from "react";
import { useTimeSeriesData } from "../../hooks/useTimeSeriesData";
import { ApexOptions } from "apexcharts";
import ReactApexChart from "react-apexcharts";
import { motion } from "framer-motion";
import { getTimeSeriesOptions } from "./chartOptions";
import { ThemeContext } from "../../contexts/ThemeContext";
import { ThemeType } from "src/types";
import { Row } from "react-bootstrap";

interface CasesTimeSeriesChartProps {
  noFilterRequest?: boolean;
  onRender: () => void;
}

const CasesTimeSeriesChart: React.FC<CasesTimeSeriesChartProps> = ({
  noFilterRequest = false,
  onRender,
}) => {
  const { theme } = React.useContext(ThemeContext);

  const { timeSeriesData, timeSeriesCount, isError } =
    useTimeSeriesData(noFilterRequest);

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
        name: "Reports",
        data: timeSeriesData,
      },
    ],
  };

  const options: ApexOptions = getTimeSeriesOptions(theme as ThemeType);

  return (
    <>
      {timeSeriesCount !== undefined && (
        <Row className={"text-center"}>
          <span className={"text-secondary"}>
            {timeSeriesCount.toLocaleString()} reports in total
          </span>
        </Row>
      )}
      <motion.div
        layout
        initial={{ opacity: 0 }}
        animate={{ opacity: 1 }}
        transition={{ duration: 0.3 }}
        exit={{ opacity: 0 }}
      >
        <ReactApexChart
          options={options}
          series={chartData.series}
          type={options.chart?.type}
        />
      </motion.div>
    </>
  );
};

export default CasesTimeSeriesChart;
