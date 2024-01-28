import CasesTimeSeriesChart from "../ApexChart/CasesTimeSeriesChart";
import TermCarousel from "../ApexChart/TermCarousel";
import "react-responsive-carousel/lib/styles/carousel.min.css";
import { Row } from "react-bootstrap";

const ChartSection = () => {
  return (
    <div className={"mt-4"}>
      <Row className={"mb-4"}>
        <CasesTimeSeriesChart />
      </Row>
      <Row className={"mb-4"}>
        <TermCarousel />
      </Row>
    </div>
  );
};

export default ChartSection;
