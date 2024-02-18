import React from "react";
import ReactApexChart from "react-apexcharts";
import { ChartDataPoint, ThemeType } from "src/types";
import { getTermCarouselOptions } from "./chartOptions";
import { ThemeContext } from "../../contexts/ThemeContext";

interface SimpleTermChartProps {
  data: ChartDataPoint[];
}

const SimpleTermChart: React.FC<SimpleTermChartProps> = ({ data }) => {
  const { theme } = React.useContext(ThemeContext);

  const chartData = {
    series: [
      {
        name: "Reports",
        data: data,
      },
    ],
  };

  const simpleChartOptions = getTermCarouselOptions(theme as ThemeType);

  return (
    <>
      <div className={"mb-3"}>
        <ReactApexChart
          options={simpleChartOptions}
          type="bar"
          series={chartData.series}
        />
      </div>
    </>
  );
};

export default SimpleTermChart;
