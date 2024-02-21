import React from "react";
import ReactApexChart from "react-apexcharts";
import { DemographicData, ThemeType } from "src/types";
import { getTermCarouselOptions } from "./chartOptions";
import { ThemeContext } from "../../contexts/ThemeContext";
import _ from "lodash";

const getChartSeries = (data: DemographicData[], aggregateType: string) => {
  const seriesDimension = aggregateType === "Age" ? "Sex" : "Age";

  return data.map((entry) => {
    return {
      name: entry.params[seriesDimension],
      data: entry.data,
    };
  }
  );
}

interface SimpleTermChartProps {
  data: DemographicData[];
  aggregateType: string;
}

const SimpleTermChart: React.FC<SimpleTermChartProps> = ({ data, aggregateType }) => {
  const { theme } = React.useContext(ThemeContext);

  const chartSeries = getChartSeries(data, aggregateType);

  const simpleChartOptions = getTermCarouselOptions(theme as ThemeType);

  const specialOptions = {
    colors: ["#59768A", "#035363", "#32B2BF", "#D5E0BE", "#CE9062", "#E0AB86", "#C7CE8A", "#6EB585", "#325951", "#6F9F9D"],
    legend: {
      show: true,
      position: "top",
    },
    chart: {
      stacked: true,
    },
    plotOptions: {
      bar: {
        horizontal: true,
        distributed: false,
      }
    },
  };

  const chartOptions = _.merge(simpleChartOptions, specialOptions);

  return (
    <>
      <div className={"mb-3"}>
        <ReactApexChart
          options={chartOptions}
          type="bar"
          series={chartSeries}
        />
      </div>
    </>
  );
};

export default SimpleTermChart;
