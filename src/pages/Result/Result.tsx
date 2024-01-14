import {useUrlParams} from "../../hooks/useUrlParams";
import { useNavigate } from "react-router-dom";
import SearchBox from "../../components/SearchBox/SearchBox";
import Header from "../../components/Header/Header";
import {Col, Container, Row} from "react-bootstrap";
import useDrugInfo from "../../hooks/useDrugInfo";
import React from "react";
import DrugPropertyBox from "../../components/DrugPropertyBox/DrugPropertyBox";

const Result = () => {
    const { params, error } = useUrlParams()
    const capitalizedTerms = params.terms.map(term => term.charAt(0).toUpperCase() + term.slice(1))
    const navigate = useNavigate()

    const drugInfo = useDrugInfo(params)

    React.useEffect(() => {
        if (error) {
            navigate('/error')
        }
    }, [error, navigate])

    React.useEffect(() => {
        if (params.terms) {
            document.title = capitalizedTerms.join(' & ') + ' - Drug Search'
        }
    }, [capitalizedTerms, params.terms])

    return (
        <>
            <Row>
                <Col className="text-center mb-5 mt-3">
                    <Row className="w-100 mx-auto">
                        <Col xs="auto" className="text-center mb-4 mt-3 mx-auto">
                            <Container className="mt-3 w-auto">
                                <Header />
                            </Container>
                        </Col>
                    </Row>
                    <SearchBox />
                </Col>
            </Row>
            <Row className="w-100 mx-auto justify-content-center">
                {drugInfo.map((drug, index) => (
                    <Col xs={12} md={6} lg={4} key={index} className="mb-4">
                        <DrugPropertyBox drug={drug} />
                    </Col>
                ))}
            </Row>
        </>

    );
}

export default Result;