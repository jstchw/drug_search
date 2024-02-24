import { Nav, Row, Col } from "react-bootstrap";
import SearchBox from "../SearchBox/SearchBox";
import Header from "../Header/Header";
import { motion } from "framer-motion";

const Scrollbar = () => {
  return (
    <motion.div
        initial={{ y: -100 }}
        animate={{ y: 0 }}
        exit={{ y: -100 }}
        transition={{ duration: 0.2 }}
        className="fixed-top py-3 m-3 rounded shadow-sm bg-body"
        >
    <Nav>
    <Row className="w-100 align-items-center">
      <Col className={"ms-4"}>
        <Header />
      </Col>

      <Col xs={4} className="d-flex justify-content-center">
        <SearchBox />
      </Col>

      <Col xs={4}></Col>
    </Row>
  </Nav>
    </motion.div>
  );
};

export default Scrollbar;
