import React from "react";
import { useTimeSeriesData } from "../../hooks/useTimeSeriesData";
import { ApexOptions } from "apexcharts";
import ReactApexChart from "react-apexcharts";
import { motion } from "framer-motion";
import { getTimeSeriesOptions } from "./chartOptions";
import { ThemeContext } from "../../contexts/ThemeContext";
import { ThemeType } from "src/types";

interface CasesTimeSeriesChartProps {
  noFilterRequest?: boolean;
  onRender: () => void;
}

const CasesTimeSeriesChart: React.FC<CasesTimeSeriesChartProps> = ({
  noFilterRequest = false,
  onRender,
}) => {
  const { theme } = React.useContext(ThemeContext);

  const { data, error } = useTimeSeriesData(noFilterRequest);

  React.useEffect(() => {
    if (data && !error) {
      onRender();
    }
  }, [data, error]);

  if (!data || error) {
    return null;
  }

  const chartData = {
    series: [
      {
        name: "Reports",
        data: data,
      },
    ],
  };

  const options: ApexOptions = getTimeSeriesOptions(theme as ThemeType);

  return (
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
  );
};

export default CasesTimeSeriesChart;
